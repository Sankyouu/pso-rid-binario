function [best_sol, best_cost, cost_history] = pso_rid_generico(funcao_custo, config_vars)

% Define a função — recebe: funcao_custo (handle da função objetivo) e config_vars (struct com tipo/min/max de cada variável)
% Retorna: best_sol (melhor solução encontrada), best_cost (melhor custo), cost_history (histórico do custo por iteração)

% SOLVER PSO RID (Real-Inteiro-Discreto) Genérico

% --- VALIDAÇÃO DE ENTRADA ---
assert(~isempty(config_vars), 'config_vars não pode ser vazio.');
% Garante que config_vars não chegou vazio — se estiver, lança erro com a mensagem indicada antes de qualquer cálculo

for i = 1:length(config_vars)        % Percorre cada variável declarada pelo usuário para verificar se os campos obrigatórios existem
    v = config_vars(i);              % Atalho local: v aponta para a i-ésima variável da config
    assert(isfield(v, 'tipo'), 'config_vars(%d) deve ter campo "tipo".', i);
    % Verifica se o campo 'tipo' existe nesta variável — sem ele o solver não sabe como tratar a dimensão

    if strcmp(v.tipo, 'R') || strcmp(v.tipo, 'I')
        % Se a variável é Real ou Inteira, precisa obrigatoriamente de 'min' e 'max'
        assert(isfield(v,'min') && isfield(v,'max'), ...
            'config_vars(%d) tipo %s requer campos "min" e "max".', i, v.tipo);
        % Lança erro claro dizendo qual variável (i) e qual tipo está com problema
    elseif strcmp(v.tipo, 'D')
        % Se é Discreta, precisa do campo 'opcoes' com pelo menos um valor
        assert(isfield(v,'opcoes') && ~isempty(v.opcoes), ...
            'config_vars(%d) tipo D requer campo "opcoes" não vazio.', i);
    else
        error('Tipo desconhecido "%s" em config_vars(%d). Use R, I ou D.', v.tipo, i);
        % Se o tipo não é R, I nem D, lança erro imediatamente informando o valor inválido recebido
    end
end

% --- 1. PARÂMETROS ---
n_particulas   = 100;     % Número de partículas no enxame — cada uma é uma solução candidata explorando o espaço
max_iter       = 1000;   % Número máximo de iterações do loop principal — limita o tempo de execução
w              = 0.8;    % Fator de inércia — controla quanto da velocidade anterior a partícula mantém (>1 explora mais, <1 converge mais)
c1             = 2.0;    % Coeficiente cognitivo — peso da atração da partícula em direção à sua própria melhor posição já vista
c2             = 2.0;    % Coeficiente social — peso da atração da partícula em direção à melhor posição já vista por todo o enxame
pm             = 0.15;   % Probabilidade de mutação — a cada iteração, 15% de chance de uma dimensão binária receber um valor aleatório
tol_estagnacao = 200;    % Critério de parada antecipada — se o GBest não melhorar por 200 iterações seguidas, encerra o loop

% --- 2. MAPEAMENTO DA PARTÍCULA (Seção 4.1, Eq. 5) ---
% Pré-alocação da struct para evitar crescimento dinâmico, decodificação binária vetorizada com produto interno
template.idx_original = 0; template.tipo = ''; template.n_bits = 0;
template.inicio = 0;       template.fim   = 0;
% Cria uma struct-modelo com todos os campos zerados — necessário para que o repmat abaixo funcione sem erro de campo heterogêneo

mapa_dimensoes = repmat(template, length(config_vars), 1);
% Pré-aloca um array de structs com o mesmo número de entradas que config_vars — evita realocação de memória a cada iteração do loop

total_dim = 0;           % Contador do tamanho total do vetor de posição de cada partícula (cresce conforme variáveis são mapeadas)

for i = 1:length(config_vars)        % Percorre cada variável declarada para construir o mapa de dimensões
    var = config_vars(i);            % Atalho local para a i-ésima variável
    info.idx_original = i;           % Guarda o índice original da variável — usado depois para acessar config_vars corretamente [FIX-2]
    info.tipo         = var.tipo;    % Copia o tipo ('R', 'I' ou 'D') para o mapa

    if strcmp(var.tipo, 'R')         % Variável Real ocupa exatamente 1 posição no vetor (valor contínuo direto)
        info.n_bits = 0;             % Reais não usam representação binária — n_bits = 0 serve como marcador
        info.inicio = total_dim + 1; % Esta dimensão começa na próxima posição disponível do vetor
        total_dim   = total_dim + 1; % Avança o contador em 1 (real ocupa 1 dimensão)
        info.fim    = total_dim;     % Esta dimensão termina na mesma posição em que começou (só 1 posição)
    else
        if strcmp(var.tipo, 'I')
            limite = var.max;        % Para inteiros, o limite é o valor máximo — precisamos de bits suficientes para representá-lo
        else  % 'D'
            limite = length(var.opcoes);   % Para discretos, o limite é a quantidade de opções disponíveis
        end
        n_bits      = max(ceil(log2(limite + 0.1)), 1);
        % Calcula quantos bits são necessários: log2(limite) arredondado para cima
        % O +0.1 evita erro de arredondamento quando limite é potência exata de 2
        % O max(..., 1) garante pelo menos 1 bit mesmo para limites muito pequenos

        info.n_bits = n_bits;        % Armazena o número de bits desta variável no mapa
        info.inicio = total_dim + 1; % Esta variável começa na próxima posição livre do vetor
        total_dim   = total_dim + n_bits;   % Avança o contador pelo número de bits desta variável
        info.fim    = total_dim;     % Esta variável termina após todos os seus bits
    end
    mapa_dimensoes(i) = info;        % Salva as informações desta variável na posição i do mapa pré-alocado
end

% --- PRÉ-COMPUTAÇÃO DE TIPOS POR DIMENSÃO ---
% Elimina o loop de busca O(n) dentro do loop principal
tipo_por_dim = cell(total_dim, 1);   % Cell array com uma entrada por dimensão do vetor — vai guardar o tipo ('R', 'I' ou 'D') de cada posição
orig_por_dim = zeros(total_dim, 1);  % Vetor numérico com uma entrada por dimensão — vai guardar o índice da variável original de cada posição

for k = 1:length(mapa_dimensoes)    % Percorre cada variável do mapa
    for d = mapa_dimensoes(k).inicio : mapa_dimensoes(k).fim
        % Percorre cada dimensão (posição do vetor) que pertence a esta variável
        tipo_por_dim{d} = mapa_dimensoes(k).tipo;        % Marca o tipo desta posição com o tipo da variável dona dela
        orig_por_dim(d) = mapa_dimensoes(k).idx_original; % Marca o índice original da variável dona desta posição
    end
end
% Resultado: dado qualquer posição d do vetor, tipo_por_dim{d} e orig_por_dim(d) respondem em O(1) sem loop de busca

% --- 3. INICIALIZAÇÃO (Seção 4.1) ---
pos = zeros(n_particulas, total_dim);   % Matriz de posições: cada linha é uma partícula, cada coluna é uma dimensão do vetor híbrido
vel = zeros(n_particulas, total_dim);   % Matriz de velocidades: mesma estrutura — começa tudo em zero

for i = 1:n_particulas               
    for k = 1:length(mapa_dimensoes) 
        info = mapa_dimensoes(k);    % Recupera as informações desta variável (tipo, início, fim, etc.)
        orig = info.idx_original;    % Índice correto da variável original em config_vars

        if strcmp(info.tipo, 'R')    % Variável Real: inicializa com valor contínuo aleatório dentro do intervalo
            vmin = config_vars(orig).min;   % Limite inferior da variável real
            vmax = config_vars(orig).max;   % Limite superior da variável real
            pos(i, info.inicio) = vmin + rand * (vmax - vmin);
            % rand gera número em [0,1]; multiplicar por (vmax-vmin) e somar vmin escala para [vmin, vmax]
        else
            pos(i, info.inicio:info.fim) = randi([0, 1], 1, info.n_bits);
            % Preenche os n_bits desta variável com 0s e 1s aleatórios — cada bit tem 50% de chance de ser 0 ou 1
        end
    end
end

% Memórias
pbest_pos   = pos;                   % PBest (Personal Best): inicializa com as posições atuais — cada partícula começa como sua própria melhor posição
pbest_custo = inf(n_particulas, 1);  % Melhor custo de cada partícula — começa em infinito pois nenhuma foi avaliada ainda
pbest_viol  = inf(n_particulas, 1);  % Melhor violação de cada partícula — começa em infinito (Método de Deb)

gbest_pos   = zeros(1, total_dim);   % GBest (Global Best): posição da melhor partícula de todo o enxame — inicializado com zeros
gbest_custo = inf;                   % Melhor custo global — começa em infinito
gbest_viol  = inf;                   % Melhor violação global — começa em infinito

cost_history    = nan(max_iter, 1);  % Vetor para gravar o custo do GBest a cada iteração — NaN indica posições ainda não preenchidas
iter_sem_melhora = 0;                % Contador de iterações consecutivas sem melhoria no GBest — usado para parada antecipada
melhor_custo_ant = inf;              % Guarda o melhor custo da iteração anterior para comparar se houve melhoria

% --- 4. LOOP PRINCIPAL ---
for iter = 1:max_iter                

    % --- A. AVALIAÇÃO (Com Método de Deb) ---
    for i = 1:n_particulas           

        vars_decodificadas = decodificar(pos(i,:), mapa_dimensoes, config_vars);
        % Converte o vetor híbrido binário/real da partícula i em valores interpretáveis pela função de custo

        [custo_atual, viol_atual] = funcao_custo(vars_decodificadas);
        % Chama a função objetivo — retorna o custo (J) e a violação das restrições desta solução

        % Atualiza PBest (3 regras de Deb)
        melhorou = false;            % Flag: assume que não houve melhoria até provar o contrário
        if     (pbest_viol(i) > 0) && (viol_atual < pbest_viol(i)), melhorou = true;
        % Regra 1: ambas as soluções são inviáveis, mas a nova viola menos — a nova é melhor
        elseif (pbest_viol(i) > 0) && (viol_atual == 0),             melhorou = true;
        % Regra 2: o PBest anterior era inviável e a nova solução é viável — a nova sempre vence
        elseif (pbest_viol(i) == 0) && (viol_atual == 0) && (custo_atual < pbest_custo(i))
            melhorou = true;
        % Regra 3: ambas são viáveis — vence a com menor custo objetivo
        end

        if melhorou                  % Se alguma das 3 regras indicou melhoria...
            pbest_custo(i)  = custo_atual;   % ...atualiza o melhor custo pessoal desta partícula
            pbest_viol(i)   = viol_atual;    % ...atualiza a melhor violação pessoal
            pbest_pos(i,:)  = pos(i,:);      % ...atualiza a melhor posição pessoal (vetor híbrido completo)
        end

        % Atualiza GBest
        melhorou_global = false;     % Flag: assume que o GBest não foi superado
        if     (gbest_viol > 0) && (pbest_viol(i) < gbest_viol), melhorou_global = true;
        % Regra 1 global: GBest inviável e o PBest desta partícula viola menos
        elseif (gbest_viol > 0) && (pbest_viol(i) == 0),          melhorou_global = true;
        % Regra 2 global: GBest inviável e o PBest desta partícula é viável
        elseif (gbest_viol == 0) && (pbest_viol(i) == 0) && (pbest_custo(i) < gbest_custo)
            melhorou_global = true;
        % Regra 3 global: ambos viáveis e o PBest desta partícula tem custo menor que o GBest atual
        end

        if melhorou_global           % Se o PBest desta partícula superou o GBest atual...
            gbest_custo = pbest_custo(i);    % ...atualiza o melhor custo global
            gbest_viol  = pbest_viol(i);     % ...atualiza a melhor violação global
            gbest_pos   = pbest_pos(i,:);    % ...atualiza a melhor posição global (o novo líder do enxame)
        end
    end

    % --- B. BUSCA LOCAL NO GBEST (Seção 5) ---    
    gbest_temp = gbest_pos;          % Copia o GBest atual para uma variável temporária — não modifica o original até confirmar melhoria
    idx_mut    = randi(total_dim);   % Sorteia aleatoriamente qual dimensão do vetor será perturbada
    tipo_mut   = tipo_por_dim{idx_mut};      % Consulta o tipo da dimensão sorteada em O(1)
    orig_mut   = orig_por_dim(idx_mut);      % Consulta o índice da variável original desta dimensão em O(1)

    if strcmp(tipo_mut, 'R')         % Se a dimensão sorteada é de uma variável Real...
        % Perturbação gaussiana (5% do intervalo)
        sigma = (config_vars(orig_mut).max - config_vars(orig_mut).min) * 0.05;
        % Desvio padrão da perturbação = 5% do intervalo total da variável — perturbação proporcional à escala do problema
        gbest_temp(idx_mut) = gbest_temp(idx_mut) + randn * sigma;
        % Adiciona ruído gaussiano à dimensão sorteada — randn() gera valor com distribuição normal (média 0, desvio 1)
        gbest_temp(idx_mut) = max(config_vars(orig_mut).min, ...
                              min(config_vars(orig_mut).max, gbest_temp(idx_mut)));
        % Clamp: garante que o valor perturbado não saia do intervalo [min, max] da variável
    else
        % Inversão de bit para dimensões binárias (I ou D)
        gbest_temp(idx_mut) = 1 - gbest_temp(idx_mut);
        % Se o bit era 0 vira 1, se era 1 vira 0 — exploração local do espaço binário
    end

    vars_mut    = decodificar(gbest_temp, mapa_dimensoes, config_vars);
    % Decodifica o GBest perturbado para avaliar se a mutação melhorou a solução
    [c_mut, v_mut] = funcao_custo(vars_mut);
    % Avalia o custo e a violação da solução perturbada

    if (v_mut == 0 && gbest_viol > 0) || (v_mut == 0 && c_mut < gbest_custo)
        % Aceita a perturbação se: o GBest antigo era inviável e o novo é viável,
        % OU ambos são viáveis mas o novo tem custo menor
        gbest_pos   = gbest_temp;    % Substitui o GBest pela versão perturbada
        gbest_custo = c_mut;         % Atualiza o custo do GBest
        gbest_viol  = v_mut;         % Atualiza a violação do GBest
    end

    cost_history(iter) = gbest_custo;   % Registra o custo do GBest desta iteração no histórico (para o gráfico de convergência)

    % Critério de parada antecipada por estagnação
    if gbest_custo < melhor_custo_ant - 1e-12
        % Se o GBest melhorou pelo menos 1e-12 em relação à iteração anterior (margem para evitar falso positivo por ruído numérico)
        melhor_custo_ant  = gbest_custo;   % Atualiza o melhor custo anterior com o novo valor
        iter_sem_melhora  = 0;             % Zera o contador de estagnação
    else
        iter_sem_melhora = iter_sem_melhora + 1;   % Incrementa o contador — mais uma iteração sem melhoria real
    end

    if mod(iter, 100) == 0 || iter == 1
        % Imprime progresso a cada 100 iterações e também na primeira — evita spam no terminal
        fprintf('Iter %4d | Violação: %.4f | Custo: %.10e | Sem melhora: %d\n', ...
                iter, gbest_viol, gbest_custo, iter_sem_melhora);
    end

    if iter_sem_melhora >= tol_estagnacao
        % Se o contador de estagnação atingiu o limite (200 iters sem melhora), encerra antecipadamente
        fprintf('>>> Parada antecipada na iteração %d (estagnação de %d iters)\n', ...
                iter, tol_estagnacao);
        cost_history = cost_history(1:iter);   % Trunca o histórico para conter só as iterações que realmente rodaram (remove NaNs do final)
        break;                       % Sai do loop principal
    end

    % --- C. ATUALIZAÇÃO DE VELOCIDADE E POSIÇÃO ---
    for i = 1:n_particulas           
        for d = 1:total_dim          

            tipo_d = tipo_por_dim{d};   % Recupera o tipo desta dimensão em O(1) — sem loop de busca
            orig_d = orig_por_dim(d);   % Recupera o índice da variável original desta dimensão em O(1)

            r1 = rand; r2 = rand;    % Dois números aleatórios independentes em [0,1] — introduzem estocasticidade na atualização
            tendencia = w * vel(i,d) ...
                + c1*r1*(pbest_pos(i,d) - pos(i,d)) ...
                + c2*r2*(gbest_pos(d)   - pos(i,d));
            % Calcula a "tendência de movimento" desta dimensão (Eq. 1 do artigo):
            %   w * vel(i,d)              → componente inercial: mantém direção anterior
            %   c1*r1*(pbest - pos)       → componente cognitiva: puxa em direção ao melhor pessoal
            %   c2*r2*(gbest - pos)       → componente social: puxa em direção ao melhor global

            if strcmp(tipo_d, 'R')   % Para dimensão de uma variável Real
                % --- Lógica Real (Seção 4.2) --- 
                vmax_d   = (config_vars(orig_d).max - config_vars(orig_d).min) / 2;
                % Velocidade máxima permitida = metade do intervalo da variável — evita que a partícula "exploda" para fora do espaço
                vel(i,d) = max(-vmax_d, min(vmax_d, tendencia));
                % Aplica clamping na velocidade: limita entre -vmax_d e +vmax_d
                pos(i,d) = pos(i,d) + vel(i,d);
                % Atualiza a posição somando a velocidade (Eq. 2 do paper)
                pos(i,d) = max(config_vars(orig_d).min, ...
                           min(config_vars(orig_d).max, pos(i,d)));
                % Clamping de posição: garante que a partícula não saia do intervalo [min, max]
            else
                % --- Lógica Binária (Seção 4.3, Eq. 8) ---
                if     tendencia > 0, v_next =  1;   % Tendência positiva → partícula deve aumentar o bit (mover para 1)
                elseif tendencia < 0, v_next = -1;   % Tendência negativa → partícula deve diminuir o bit (mover para 0)
                else,                 v_next =  0;   % Tendência nula → partícula permanece no mesmo bit
                end

                % Mutação estocástica
                if rand < pm         % Com probabilidade pm (15%), aplica mutação — introduz diversidade para escapar de mínimos locais
                    v_next = [-1, 0, 1];     % Define os 3 movimentos possíveis
                    v_next = v_next(randi(3));   % Escolhe um dos 3 aleatoriamente com igual probabilidade
                end

                vel(i,d) = v_next;           % Atualiza a velocidade binária: -1, 0 ou +1
                pos(i,d) = pos(i,d) + v_next;   % Atualiza a posição somando o movimento binário

                % Clamp binário
                pos(i,d) = max(0, min(1, pos(i,d)));
                % Garante que o bit permaneça em {0, 1}: se passou de 1 vira 1, se ficou abaixo de 0 vira 0
            end
        end
    end

end  % fim do loop principal

best_cost = gbest_custo;             % Copia o melhor custo global encontrado para a variável de saída
best_sol  = decodificar(gbest_pos, mapa_dimensoes, config_vars);
% Decodifica a melhor posição global para o espaço de variáveis originais e retorna como solução final

end  % fim da função principal


% DECODIFICADOR (Eq. 6)
% Converte o vetor híbrido em variáveis interpretáveis

function vars = decodificar(vetor_hibrido, mapa, config)
% Recebe o vetor híbrido (posição bruta da partícula), o mapa de dimensões e a config original
% Retorna vars: vetor com os valores reais/inteiros/discretos de cada variável

n_vars = length(config);             % Número total de variáveis do problema
vars   = zeros(1, n_vars);           % Pré-aloca o vetor de saída com zeros

for k = 1:n_vars                     
    info = mapa(k);                  % Recupera as informações de mapeamento desta variável (tipo, início, fim)

    if strcmp(info.tipo, 'R')        % Variável Real: o valor está diretamente na posição do vetor, sem decodificação
        vars(k) = vetor_hibrido(info.inicio);   % Copia o valor contínuo direto da posição correspondente

    else                             % Variável Inteira ou Discreta: precisa converter bits para número
        bits = vetor_hibrido(info.inicio : info.fim);   % Extrai os bits desta variável do vetor híbrido
        n    = length(bits);         % Número de bits desta variável

        % Conversão vetorizada: produto interno com pesos binários
        pesos   = 2 .^ ((n-1):-1:0);   % Gera os pesos posicionais: [2^(n-1), 2^(n-2), ..., 2, 1]
        % Exemplo para 3 bits: pesos = [4, 2, 1] — bit mais significativo tem peso 4
        decimal = round(bits) * pesos';
        % round(bits): arredonda cada bit para 0 ou 1 (necessário pois valores como 0.99 surgem durante a atualização)
        % * pesos': produto interno — soma cada bit multiplicado pelo seu peso → número decimal

        if strcmp(info.tipo, 'I')    % Variável Inteira: o decimal é o valor diretamente
            vars(k) = min(decimal, config(k).max);
            % Clamp no máximo: se os bits representarem um valor maior que o permitido, usa o máximo

        else  % 'D'                  % Variável Discreta: o decimal é um índice na lista de opções
            % mod() garante que índice 0 seja remapeado para 1
            N       = length(config(k).opcoes);   % Número de opções disponíveis para esta variável
            idx     = mod(decimal, N) + 1;        % mod mapeia qualquer decimal para [0, N-1], o +1 desloca para [1, N] — sempre válido em MATLAB
            vars(k) = config(k).opcoes(idx);      % Seleciona o valor correspondente ao índice calculado
        end
    end
end
end