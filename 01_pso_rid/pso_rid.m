function [best_sol, best_cost, cost_history, details] = pso_rid(funcao_custo, config_vars, pso_params)
% PSO_RID  Particle Swarm Optimization Real-Inteiro-Discreto (codificacao binaria).
%
% -------------------------------------------------------------------------
% MAPA DESTE ARQUIVO  (para navegar, busque pela marca ">>>")
% -------------------------------------------------------------------------
%
%   >>> ACESSO AOS AUXILIARES (teste)
%   >>> VALIDACAO DE ENTRADA
%   >>> HIPERPARAMETROS
%   >>> MAPEAMENTO DA PARTICULA .................. Eq. (5)
%   >>> INICIALIZACAO DO ENXAME
%   >>> LOOP PRINCIPAL
%         A. AVALIACAO E ATUALIZACAO DE p-BEST / g-BEST ... Eq. (3)-(4)
%         B. BUSCA LOCAL NO g-BEST ....................... Sec. 5
%         C. ATUALIZACAO DE VELOCIDADE E POSICAO ......... Eq. (1)-(2)
%   >>> SAIDAS
%   >>> LOCAL mapear_dimensoes ....... Eq. (5)      dimensoes da particula
%   >>> LOCAL decodificar ............ Eq. (6)      binario -> real/int/discreto
%   >>> LOCAL velocidade_binaria ..... Eq. (7)-(8)  velocidade binaria
%   >>> LOCAL mutacao_polinomial ..... Eq. (9)      mutacao polinomial
%   >>> LOCAL domina_deb ............. pag. 316     regra de viabilidade de Deb
%   >>> LOCAL amostrar_particula ......... Sec. 4.1     inicializacao aleatoria
%   >>> LOCAL get_param .................. utilitario   default de hiperparametro
%
% -------------------------------------------------------------------------
% REFERENCIAS TEORICAS
% -------------------------------------------------------------------------
%
% [DF2011] Datta, D.; Figueira, J.R.
%          "A real-integer-discrete-coded particle swarm optimization for
%           design problems". Applied Soft Computing 11 (2011) 3625-3633.
%          -> 00_docs/artigos/Areal-integer-discrete-coded particle swarm
%             optimization for design problems.pdf
%
% [DEB2000] Deb, K. "An efficient constraint handling method for genetic
%           algorithms". Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%          -> 00_docs/artigos/An efficient constraint handling method for
%             genetic algorithms.pdf
%
% Mapeamento equacao:
%   Eq. (1) [DF2011] velocidade real ............... secao C
%   Eq. (2) [DF2011] posicao real .................. secao C
%   Eq. (3) [DF2011] atualizacao do p-best ......... secao A
%   Eq. (4) [DF2011] atualizacao do g-best ......... secao A
%   Eq. (5) [DF2011] dimensoes da particula ........ mapear_dimensoes
%   Eq. (6) [DF2011] decodificacao binaria ......... decodificar
%   Eq. (7)-(8), Tab. 1-2 [DF2011] veloc. binaria .. velocidade_binaria
%   Eq. (9) [DF2011] mutacao polinomial ............ mutacao_polinomial
%   Sec. 5  [DF2011] busca local no g-best ......... secao B
%   pag. 316 [DEB2000] dominancia (3 criterios) .... domina_deb
%
% -------------------------------------------------------------------------
% ENTRADAS
% -------------------------------------------------------------------------
%
%   funcao_custo : handle @(x) [custo, violacao]
%   config_vars  : struct array definindo cada variavel:
%                    .tipo = 'R' (real)     -> requer .min e .max
%                    .tipo = 'I' (inteiro)  -> requer .min e .max
%                    .tipo = 'D' (discreto) -> requer .opcoes (vetor)
%   pso_params   : (opcional) struct de hiperparametros
%
% SAIDAS
%   best_sol     : melhor solucao decodificada
%   best_cost    : custo objetivo da melhor solucao
%   cost_history : historico de convergencia do g-best (vetor coluna).
%                  ATENCAO: NAO e monotonico. Pelo criterio 1 de [DEB2000]
%                  um projeto viavel e preferido a qualquer inviavel, entao
%                  o custo SOBE no instante em que o enxame troca um
%                  inviavel-barato por um viavel-caro. O historico que e
%                  monotonico nao crescente e details.viol_history.
%   details      : struct diagnostico (inclui .viol_history)
%
% -------------------------------------------------------------------------
% ACESSO AOS AUXILIARES PARA TESTE
% -------------------------------------------------------------------------
%
% A chamada especial abaixo devolve os handles para os testes unitários:
%
%   aux = pso_rid('auxiliares');
%   aux.mapear_dimensoes    -> mapear_dimensoes
%   aux.decodificar         -> decodificar
%   aux.velocidade_binaria  -> velocidade_binaria
%   aux.mutacao_polinomial  -> mutacao_polinomial
%   aux.domina_deb          -> domina_deb
%
% -------------------------------------------------------------------------
% >>> ACESSO AOS AUXILIARES (teste)
% -------------------------------------------------------------------------

if nargin == 1 && ischar(funcao_custo) && strcmp(funcao_custo, 'auxiliares')
    best_sol = struct( ...
        'mapear_dimensoes',   @mapear_dimensoes,   ...
        'decodificar',        @decodificar,        ...
        'velocidade_binaria', @velocidade_binaria, ...
        'mutacao_polinomial', @mutacao_polinomial, ...
        'domina_deb',         @domina_deb);
    best_cost    = [];
    cost_history = [];
    details      = struct();
    return;
end

% -------------------------------------------------------------------------
% >>> VALIDACAO DE ENTRADA
% -------------------------------------------------------------------------

assert(~isempty(config_vars), 'config_vars nao pode ser vazio.');

for i = 1:length(config_vars)
    v = config_vars(i);
    assert(isfield(v, 'tipo'), 'config_vars(%d) deve ter campo "tipo".', i);

    if strcmp(v.tipo, 'R') || strcmp(v.tipo, 'I')
        assert(isfield(v,'min') && isfield(v,'max'), ...
            'config_vars(%d) tipo %s requer campos "min" e "max".', i, v.tipo);
        assert(v.max > v.min, ...
            'config_vars(%d): "max" deve ser maior que "min".', i);
    elseif strcmp(v.tipo, 'D')
        assert(isfield(v,'opcoes') && ~isempty(v.opcoes), ...
            'config_vars(%d) tipo D requer campo "opcoes" nao vazio.', i);
    else
        error('Tipo desconhecido "%s" em config_vars(%d). Use R, I ou D.', v.tipo, i);
    end
end

if nargin < 3 || isempty(pso_params)
    pso_params = struct();
end

campos_reconhecidos = {'n_particulas', 'max_iter', 'w', 'c1', 'c2', 'pm', ...
                       'auto_adaptativo', 'eta_min', 'eta_max', 'p_busca_local', ...
                       'tol_estagnacao', 'reinit_freq', 'reinit_pct', ...
                       'decodificacao_discreta', ...
                       'verbose', 'print_interval'};
campos_desconhecidos = setdiff(fieldnames(pso_params), campos_reconhecidos);
if ~isempty(campos_desconhecidos)
    warning('pso_rid:parametroDesconhecido', ...
        ['Campo(s) de pso_params nao reconhecido(s) e IGNORADO(S): %s. ' ...
         'Campos validos: %s.'], ...
        strjoin(campos_desconhecidos, ', '), strjoin(campos_reconhecidos, ', '));
end

% -------------------------------------------------------------------------
% >>> HIPERPARAMETROS
% -------------------------------------------------------------------------
% [DF2011] Sec. 5, pag. 3628: "In each run, the mutation probability p_m for
% binary variables is assigned the initial value of 15%. For handling real
% variables, the assigned initial values of the inertia constant w, cognitive
% behavioral factor c1, and social behavioral factor c2 are 1.0, 1.0 and 2.0,
% respectively."
% -------------------------------------------------------------------------

n_particulas    = get_param(pso_params, 'n_particulas',   100);
max_iter        = get_param(pso_params, 'max_iter',      1000);
w0              = get_param(pso_params, 'w',              1.0);  % [DF2011] Sec. 5
c1_0            = get_param(pso_params, 'c1',             1.0);  % [DF2011] Sec. 5
c2_0            = get_param(pso_params, 'c2',             2.0);  % [DF2011] Sec. 5
pm0             = get_param(pso_params, 'pm',            0.15);  % [DF2011] Sec. 4.3/5

% [DF2011] Sec. 5: "each of p_m, w, c1 and c2 is made instant-wise
% self-adaptive in a range from zero to its assigned initial value."

auto_adaptativo = get_param(pso_params, 'auto_adaptativo', true);

% [DF2011] Sec. 5: indice de distribuicao da mutacao polinomial (Eq. 9),
% "assigned a random value, in the range of [25,45], in different runs"

eta_min         = get_param(pso_params, 'eta_min',         25);
eta_max         = get_param(pso_params, 'eta_max',         45);

% [DF2011] Sec. 5 (busca local): "a binary position of a particle is altered
% ... with a small random probability in the range of ]0,0.05["

p_busca_local   = get_param(pso_params, 'p_busca_local',  0.05);

% -------------------------------------------------------------------------
% ESTRATEGIA DE DECODIFICACAO DISCRETA
% -------------------------------------------------------------------------
%   'proporcional' (PADRAO) indice = floor(x * N / 2^B) + 1.
%                  Todos os 2^B codigos caem no catalogo: nao ha codigo
%                  morto. Custo: vies estatistico — com N=13 em 4 bits, as
%                  opcoes 1, 5 e 9 recebem 2 codigos cada (12.5% de chance)
%                  contra 1 codigo das outras dez (6.25%).
%
%   'datta'        indice = x + 1, fiel a [DF2011] Sec. 4.1. Indices > N nao
%                  existem no catalogo e viram inviabilidade estrutural.
%
%   Com N=13 e B=4, sobram 3 dos 16 codigos por variavel. Com 10 variaveis
%   isso compoe: (13/16)^10 = 12.5% — apenas uma particula em oito e
%   inteiramente valida, e as outras 87.5% nunca viram projeto avaliavel.
%   [DF2011] nunca testou esse regime: no seu unico exemplo com discreta
%   (Sec. 5.2, mola) ha UMA variavel, e o desperdicio nao compoe.
%
%   Evidencia (30 sementes pareadas, orcamento igualado em avaliacoes FEM,
%   benchmark Hadi linear — log em 04_resultados/logs/):
%
%     configuracao                 melhor     media    desvio
%     'datta'                     2357.97   2572.64    173.05
%     'proporcional'              2342.11   2504.72     97.40
%
%   Sem codigos mortos, 100% das particulas viram
%   projeto avaliavel contra 12.5% no 'datta'. O enxame carrega cerca de
%   oito vezes mais informacao de projeto por iteracao.

decod_discreta  = get_param(pso_params, 'decodificacao_discreta', 'proporcional');
assert(any(strcmp(decod_discreta, {'datta', 'proporcional'})), ...
    'decodificacao_discreta deve ser ''datta'' ou ''proporcional''. Recebido: ''%s''.', ...
    decod_discreta);

tol_estagnacao  = get_param(pso_params, 'tol_estagnacao', 200);
reinit_freq     = get_param(pso_params, 'reinit_freq',    100);
reinit_pct      = get_param(pso_params, 'reinit_pct',    0.10);
verbose         = get_param(pso_params, 'verbose',       true);
print_interval  = get_param(pso_params, 'print_interval', 100);

% eta sorteado uma vez por execucao ("in different runs of a problem")

eta_run = eta_min + rand * (eta_max - eta_min);

% -------------------------------------------------------------------------
% >>> MAPEAMENTO DA PARTICULA — Eq. (5) [DF2011]
% -------------------------------------------------------------------------

[mapa_dimensoes, total_dim] = mapear_dimensoes(config_vars);

% Lookup O(1) de tipo / variavel de origem por dimensao

tipo_por_dim = cell(total_dim, 1);
orig_por_dim = zeros(total_dim, 1);
for k = 1:length(mapa_dimensoes)
    for d = mapa_dimensoes(k).inicio : mapa_dimensoes(k).fim
        tipo_por_dim{d} = mapa_dimensoes(k).tipo;
        orig_por_dim(d) = mapa_dimensoes(k).idx_original;
    end
end

% -------------------------------------------------------------------------
% >>> INICIALIZACAO DO ENXAME
% -------------------------------------------------------------------------
% [DF2011] Sec. 4.1: "A real dimension is initialized by a random real value
% in the given range of the real variable which is represented by that
% dimension, while all the binary dimensions are initialized randomly by 0 or
% 1 with 50% probability. On the other hand, all the velocity components are
% assigned the initial value of 0."
% -------------------------------------------------------------------------

pos = zeros(n_particulas, total_dim);
vel = zeros(n_particulas, total_dim);   % velocidades iniciam em zero

for i = 1:n_particulas
    pos(i,:) = amostrar_particula(mapa_dimensoes, config_vars, total_dim);
end

pbest_pos   = pos;
pbest_custo = inf(n_particulas, 1);
pbest_viol  = inf(n_particulas, 1);

gbest_pos   = zeros(1, total_dim);
gbest_custo = inf;
gbest_viol  = inf;

cost_history     = nan(max_iter, 1);

% Historico da VIOLACAO do g-best. Necessario para interpretar
% cost_history: sob [DEB2000] o custo do g-best NAO e monotonico — no
% instante em que o enxame passa de inviavel-barato para viavel-caro, o
% criterio 1 manda trocar e o custo SOBE. Ja a violacao e monotonica nao
% crescente (criterios 1 e 3), e e esse o invariante verificavel.

viol_history     = nan(max_iter, 1);
iter_sem_melhora = 0;
melhor_custo_ant = inf;
iter_final       = max_iter;
n_avaliacoes     = 0;
n_overflow       = 0;

% -------------------------------------------------------------------------
% >>> LOOP PRINCIPAL
% -------------------------------------------------------------------------

for iter = 1:max_iter

    % --- Auto-adaptacao instantanea dos parametros ([DF2011] Sec. 5) ---
    if auto_adaptativo
        w_atual  = rand * w0;
        c1_atual = rand * c1_0;
        c2_atual = rand * c2_0;
        pm_atual = rand * pm0;
    else
        w_atual  = w0;
        c1_atual = c1_0;
        c2_atual = c2_0;
        pm_atual = pm0;
    end

    % ---------------------------------------------------------------------
    % A. AVALIACAO E ATUALIZACAO DE p-BEST / g-BEST
    %    Eq. (3) e (4) [DF2011], dominancia de [DEB2000] pag. 316
    % ---------------------------------------------------------------------
    
    for i = 1:n_particulas
        [vars_i, viol_estrutural] = decodificar(pos(i,:), mapa_dimensoes, ...
                                                    config_vars, decod_discreta);

        if viol_estrutural > 0

            % Valor fora do dominio (ex.: indice inexistente no catalogo):
            % a solucao nao pode ser avaliada fisicamente. Tratada como
            % inviavel, sem calcular a funcao objetivo — coerente com
            % [DEB2000] pag. 316.
            %
            % A violacao vale Inf, e nao a magnitude do estouro: uma
            % particula fora do catalogo nao e um projeto ruim, e a AUSENCIA
            % de projeto, e deve perder para qualquer projeto avaliavel. As
            % tres regras de [DEB2000] ja produzem isso sozinhas —
            % Inf < Inf e falso (nenhuma preferencia entre estouros) e
            % 1e8 < Inf e verdadeiro (avaliavel vence estouro).

            custo_atual = inf;
            viol_atual  = inf;
            n_overflow  = n_overflow + 1;
        else
            [custo_atual, viol_atual] = funcao_custo(vars_i);
            n_avaliacoes = n_avaliacoes + 1;
        end

        % Eq. (3) [DF2011]: atualizacao do p-best
        if domina_deb(custo_atual, viol_atual, pbest_custo(i), pbest_viol(i))
            pbest_custo(i) = custo_atual;
            pbest_viol(i)  = viol_atual;
            pbest_pos(i,:) = pos(i,:);
        end

        % Eq. (4) [DF2011]: atualizacao do g-best
        if domina_deb(pbest_custo(i), pbest_viol(i), gbest_custo, gbest_viol)
            gbest_custo = pbest_custo(i);
            gbest_viol  = pbest_viol(i);
            gbest_pos   = pbest_pos(i,:);
        end
    end

    % ---------------------------------------------------------------------
    % B. BUSCA LOCAL NO g-BEST
    %    [DF2011] Sec. 5: "a local search scheme is applied to the g-best
    %    particle in order to avoid any local optimum by exploring its
    %    neighborhood. In this scheme, a binary position of a particle is
    %    altered, either from 0 to 1 or from 1 to 0, with a small random
    %    probability in the range of ]0,0.05[. Similarly, a real position is
    %    also changed by applying Eq. (9)."
    %
    %    [D5] A versao anterior perturbava UMA unica dimensao sorteada, com
    %    ruido gaussiano nas reais.
    % ---------------------------------------------------------------------

    if ~isinf(gbest_custo) || ~isinf(gbest_viol)
        gbest_temp = gbest_pos;
        p_flip     = rand * p_busca_local;    % probabilidade em ]0, 0.05[

        for d = 1:total_dim
            orig_d = orig_por_dim(d);
            if strcmp(tipo_por_dim{d}, 'R')
                gbest_temp(d) = mutacao_polinomial(gbest_temp(d), ...
                                    config_vars(orig_d).min, ...
                                    config_vars(orig_d).max, eta_run);
            else
                if rand < p_flip
                    gbest_temp(d) = 1 - gbest_temp(d);
                end
            end
        end

        [vars_mut, viol_mut_estr] = decodificar(gbest_temp, mapa_dimensoes, ...
                                                    config_vars, decod_discreta);
        if viol_mut_estr > 0
            c_mut = inf;
            v_mut = inf;   % fora do catalogo — ver secao A
            n_overflow = n_overflow + 1;
        else
            [c_mut, v_mut] = funcao_custo(vars_mut);
            n_avaliacoes = n_avaliacoes + 1;
        end

        if domina_deb(c_mut, v_mut, gbest_custo, gbest_viol)
            gbest_pos   = gbest_temp;
            gbest_custo = c_mut;
            gbest_viol  = v_mut;
        end
    end

    cost_history(iter) = gbest_custo;
    viol_history(iter) = gbest_viol;

    % --- Controle de estagnacao ([DF2011] Sec. 5: "a run is terminated in
    %     between when no improvement in the best objective value is noticed")

    if gbest_custo < melhor_custo_ant - 1e-12
        melhor_custo_ant = gbest_custo;
        iter_sem_melhora = 0;
    else
        iter_sem_melhora = iter_sem_melhora + 1;
    end

    if verbose && (mod(iter, print_interval) == 0 || iter == 1)
        fprintf(['Iter %4d (w=%.3f c1=%.3f c2=%.3f pm=%.3f) | ' ...
                 'Violacao: %.4e | Custo: %.6e | Sem melhora: %d\n'], ...
                iter, w_atual, c1_atual, c2_atual, pm_atual, ...
                gbest_viol, gbest_custo, iter_sem_melhora);
    end

    % --- Reinicializacao de particulas estagnadas -------------------------
    % NOTA: mecanismo NAO descrito em [DF2011]; heuristica de diversidade
    % herdada da versao anterior. Use reinit_freq = 0 para o comportamento
    % estrito do artigo.

    if reinit_freq > 0 && mod(iter, reinit_freq) == 0 && iter_sem_melhora > 20
        n_reinit = max(1, round(reinit_pct * n_particulas));
        [~, idx_ordenados] = sort(pbest_viol, 'descend');
        for r_idx = 1:n_reinit
            p_idx = idx_ordenados(r_idx);
            pos(p_idx,:)       = amostrar_particula(mapa_dimensoes, config_vars, total_dim);
            vel(p_idx,:)       = 0;
            pbest_custo(p_idx) = inf;
            pbest_viol(p_idx)  = inf;
        end
    end

    if iter_sem_melhora >= tol_estagnacao
        if verbose
            fprintf('>>> Parada antecipada na iteracao %d (estagnacao de %d iteracoes)\n', ...
                    iter, tol_estagnacao);
        end
        cost_history = cost_history(1:iter);
        viol_history = viol_history(1:iter);
        iter_final   = iter;
        break;
    end

    % ---------------------------------------------------------------------
    % C. ATUALIZACAO DE VELOCIDADE E POSICAO
    % ---------------------------------------------------------------------

    for i = 1:n_particulas
        for d = 1:total_dim
            orig_d = orig_por_dim(d);
            r1 = rand; r2 = rand;

            % Eq. (1) [DF2011]:
            %   v = w*v + c1*r1*(pbest - x) + c2*r2*(gbest - x)
            % Para dimensoes binarias, este valor real e a "tendencia" que
            % sera mapeada para {-1,0,+1} pela Eq. (8).

            tendencia = w_atual * vel(i,d) ...
                      + c1_atual * r1 * (pbest_pos(i,d) - pos(i,d)) ...
                      + c2_atual * r2 * (gbest_pos(d)   - pos(i,d));

            if strcmp(tipo_por_dim{d}, 'R')

                % ---------- VARIAVEL REAL ----------

                vmin_d = config_vars(orig_d).min;
                vmax_d = config_vars(orig_d).max;

                % --- LIMITE DE VELOCIDADE (V_max) ---
                % Nao consta de [DF2011], mas e pratica padrao do PSO
                % continuo desde Shi & Eberhart (1998), "A modified
                % particle swarm optimizer", IEEE ICEC, p. 69-73: sem
                % teto, o termo de inercia w*v acumula e a particula
                % passa a saltar de um extremo do dominio ao outro,
                % perdendo a capacidade de refinar. A recomendacao usual
                % e V_max = k*(x_max - x_min) com k em [0.1, 1.0]; usamos
                % k = 0.5.

                v_lim    = (vmax_d - vmin_d) / 2;
                vel(i,d) = max(-v_lim, min(v_lim, tendencia));

                % Eq. (2) [DF2011]: x^(t+1) = x^(t) + v^(t+1)
                pos(i,d) = pos(i,d) + vel(i,d);

                % --- CLAMP DE POSICAO ---
                % Obrigatorio, nao opcional. Sem ele, uma dimensao real
                % fora de [min, max] quebra dois consumidores a jusante:
                %   (a) mutacao_polinomial eleva (x_u - x)/(x_u - x_l) a
                %       potencia (eta+1) com eta REAL; base negativa e
                %       expoente fracionario dao resultado COMPLEXO, que
                %       o min/max final daquela funcao nao detecta
                %       (MATLAB compara complexos por modulo);
                %   (b) decodificar repassa o valor cru a funcao de
                %       custo — uma area negativa produz matriz de
                %       rigidez sem sentido fisico no FEM.
                % [DF2011] trata limites como restricao no caso INTEIRO
                % (ver [D2]), mas ali o valor cru ainda e avaliavel; aqui
                % nao e, entao a saturacao e a unica saida coerente.

                pos(i,d) = max(vmin_d, min(vmax_d, pos(i,d)));

            else
                % ---------- VARIAVEL BINARIA ----------
                % [D3][D4] Eq. (7)-(8) e Tabelas 1-2 [DF2011]
                [v_next, x_next] = velocidade_binaria(pos(i,d), tendencia, pm_atual);
                vel(i,d) = v_next;
                pos(i,d) = x_next;
            end
        end
    end
end

% -------------------------------------------------------------------------
% >>> SAIDAS
% -------------------------------------------------------------------------

best_cost = gbest_custo;
[best_sol, viol_final_estrutural] = decodificar(gbest_pos, mapa_dimensoes, ...
                                                    config_vars, decod_discreta);

details.iter_executadas = iter_final;
details.gbest_viol      = gbest_viol;
details.gbest_raw_pos   = gbest_pos;
details.mapa_dimensoes  = mapa_dimensoes;
details.total_dim       = total_dim;
details.viol_history    = viol_history;
details.n_avaliacoes    = n_avaliacoes;
details.n_overflow      = n_overflow;
details.eta_run         = eta_run;
details.viol_estrutural = viol_final_estrutural;
details.decodificacao   = decod_discreta;

end


% #########################################################################
% #########################################################################
% ##                                                                     ##
% ##                          FUNCOES LOCAIS                             ##
% ##                                                                     ##
% ##  Uma por formula de [DF2011] / [DEB2000], na ordem das equacoes.     ##
% ##  Sao privadas a este arquivo; para exercita-las isoladamente use     ##
% ##  aux = pso_rid('auxiliares') (ver cabecalho).                        ##
% ##                                                                     ##
% #########################################################################
% #########################################################################


% -------------------------------------------------------------------------
% >>> mapear_dimensoes — Eq. (5) [DF2011]
% -------------------------------------------------------------------------

function [mapa, total_dim] = mapear_dimensoes(config_vars)

% RID_MAPEAR_DIMENSOES  Define quantas dimensoes cada variavel ocupa na particula.
%
% Eq. (5) [DF2011]:   L = R + sum_j B_j^(I) + sum_k B_k^(D)
%   L       = numero total de dimensoes posicionais da particula
%   R       = numero de variaveis reais (1 dimensao cada)
%   B_j^(I) = bits da j-esima variavel inteira
%   B_k^(D) = bits da k-esima variavel discreta
%
% [DF2011] Sec. 4.1:
%   "The number of dimensions required for an integer variable is determined
%    from the maximum number of binary bits required to represent its upper
%    limit. It is to be mentioned that a discrete variable is also dealt with
%    as an integer variable, an integer value of which represents the index of
%    its actual discrete value, so that the lower limit of such a variable is 1
%    and the upper limit is the number of its allowable discrete values."
%
% DIMENSIONAMENTO DE BITS (deduzido da Eq. 6 e conferido com os exemplos
% numericos do proprio artigo):
%
%   INTEIRO  : a Eq. (6) precisa gerar o proprio valor var.max, logo
%              2^B - 1 >= max  =>  B = ceil(log2(max + 1)).
%              Confere com Sec. 5.1 (trem de engrenagens): 12 <= z <= 60
%              com "6 binary bits" -> ceil(log2(61)) = 6.
%              Confere com Sec. 5.2 (mola): N inteiro em [1,70] com
%              "7 binary bits" -> ceil(log2(71)) = 7.
%
%   DISCRETO : o indice vai de 1 a N, ou seja x_int (Eq. 6) vai de 0 a N-1,
%              logo 2^B >= N  =>  B = ceil(log2(N)).
%              Confere com Sec. 5.2 (mola): 42 diametros discretos com
%              "6 binary bits" -> ceil(log2(42)) = 6.
%
% ENTRADA
%   config_vars : struct array com .tipo ('R','I','D') e os campos exigidos
%                 por cada tipo (.min/.max para R e I; .opcoes para D)
%
% SAIDAS
%   mapa      : struct array com, para cada variavel,
%                 .idx_original, .tipo, .n_bits, .inicio, .fim
%   total_dim : numero total de dimensoes da particula (L da Eq. 5)

template.idx_original = 0;
template.tipo         = '';
template.n_bits       = 0;
template.inicio       = 0;
template.fim          = 0;

mapa      = repmat(template, length(config_vars), 1);
total_dim = 0;

for i = 1:length(config_vars)
    var = config_vars(i);

    info.idx_original = i;
    info.tipo         = var.tipo;

    if strcmp(var.tipo, 'R')
        % Variavel real: 1 dimensao, sem bits. Eq. (5) [DF2011], termo R.
        info.n_bits = 0;
        info.inicio = total_dim + 1;
        total_dim   = total_dim + 1;
        info.fim    = total_dim;
    else
        if strcmp(var.tipo, 'I')
            info.n_bits = max(ceil(log2(var.max + 1)), 1);
        else  % 'D'
            info.n_bits = max(ceil(log2(length(var.opcoes))), 1);
        end
        info.inicio = total_dim + 1;
        total_dim   = total_dim + info.n_bits;
        info.fim    = total_dim;
    end

    mapa(i) = info;
end

end


% -------------------------------------------------------------------------
% >>> decodificar — Eq. (6) [DF2011]
% -------------------------------------------------------------------------

function [vars, viol_estrutural] = decodificar(vetor_hibrido, mapa, config, modo_discreto)
% RID_DECODIFICAR  Converte a particula hibrida (real + binaria) em valores de projeto.
%
% -------------------------------------------------------------------------
% Eq. (6) [DF2011]:   x = sum_{i=1}^{B} 2^(B-i) * b_i
%   x   = valor inteiro da variavel
%   B   = numero de bits
%   b_i = i-esimo bit
% (codificacao binaria natural, bit mais significativo a esquerda)
% -------------------------------------------------------------------------
%
% [D1] VARIAVEL DISCRETA  — decisao registrada em 2026-08-31
%   [DF2011] Sec. 4.1 trata a discreta como inteira cujo valor representa o
%   INDICE do valor discreto real, com limite inferior 1 e limite superior N.
%
%   modo_discreto = 'proporcional'  (PADRAO; NAO consta de [DF2011])
%       indice = floor(x * N / 2^B) + 1. Elimina 100% do desperdicio de
%       codigos, ao custo de um vies estatistico: com N=13 em B=4 bits, as
%       opcoes 1, 5 e 9 recebem 2 dos 16 codigos cada e as outras dez
%       recebem 1. Escolhido por ser significativamente mais consistente —
%       ver a justificativa completa e os numeros em [D1], na secao de
%       hiperparametros no inicio deste arquivo.
%
%   modo_discreto = 'datta'  (fiel ao artigo)
%       indice = x + 1, valido em 1..N. Indices > N sao ESTRUTURALMENTE
%       INVIAVEIS (o valor nao existe no catalogo) e geram viol_estrutural.
%       Consequencia: para N que nao seja potencia de 2, a fracao
%       (2^B - N)/2^B dos codigos por variavel e desperdicada, e o
%       desperdicio COMPOE entre variaveis — com 10 variaveis discretas de
%       N=13, apenas 12.5% das particulas sao inteiramente validas.
%
% [D2] VARIAVEL INTEIRA
%   O valor decodificado e usado DIRETAMENTE (Eq. 6), sem saturacao interna.
%   Em [DF2011] os limites do inteiro sao tratados como RESTRICAO do problema
%   (ver Eq. 12, trem de engrenagens: 12 <= z <= 60 com 6 bits => 0..63), e
%   nao por clamp no decodificador.
%   A versao anterior deste solver fazia min(x, max) e ignorava completamente
%   o limite inferior .min — bug corrigido aqui.
%
% -------------------------------------------------------------------------

% TRATAMENTO DE VALORES FORA DO DOMINIO
%   Sao sinalizados via viol_estrutural > 0. O chamador deve trata-los pela
%   regra de [DEB2000] (comparacao apenas por violacao), SEM avaliar a funcao
%   objetivo — coerente com [DEB2000] pag. 316: "It does not make sense to
%   compute the objective function value of an infeasible solution".
%
%   QUANDO ISSO OCORRE. Apenas por dois caminhos:
%     - variavel 'I' (inteira) fora de [min, max]  — sempre ativo;
%     - variavel 'D' (discreta) em modo 'datta'    — indice > N.
%   No modo 'proporcional' (o PADRAO) o indice e saturado em 1..N por
%   construcao e o ramo discreto NUNCA contribui.

% -------------------------------------------------------------------------
%
% ENTRADAS
%   vetor_hibrido : vetor 1 x total_dim com a posicao da particula
%   mapa          : saida de mapear_dimensoes
%   config        : struct array de configuracao das variaveis
%   modo_discreto : (opcional) 'proporcional' (padrao) ou 'datta'
%
% SAIDAS
%   vars            : vetor 1 x n_vars com os valores de projeto decodificados
%   viol_estrutural : soma das violacoes de dominio (0 = solucao avaliavel)

if nargin < 4 || isempty(modo_discreto)
    modo_discreto = 'proporcional';   % ver [D1]
end

n_vars          = length(config);
vars            = zeros(1, n_vars);
viol_estrutural = 0;

for k = 1:n_vars
    info = mapa(k);

    if strcmp(info.tipo, 'R')
        % Variavel real: valor lido diretamente da dimensao real
        vars(k) = vetor_hibrido(info.inicio);

    else
        bits = round(vetor_hibrido(info.inicio : info.fim));
        n    = length(bits);

        % --- Eq. (6) [DF2011]: pesos 2^(B-i), i = 1..B ---
        pesos = 2 .^ ((n-1):-1:0);
        x_int = bits(:)' * pesos(:);

        if strcmp(info.tipo, 'I')
            % ---------- [D2] INTEIRO ----------
            vars(k) = x_int;
            if x_int < config(k).min
                viol_estrutural = viol_estrutural + (config(k).min - x_int);
            elseif x_int > config(k).max
                viol_estrutural = viol_estrutural + (x_int - config(k).max);
            end

        else
            % ---------- [D1] DISCRETO ----------
            N = length(config(k).opcoes);

            if strcmp(modo_discreto, 'datta')
                % Fiel a [DF2011] Sec. 4.1: indice = x + 1, limites 1..N
                idx = x_int + 1;
                if idx > N
                    viol_estrutural = viol_estrutural + (idx - N);
                    vars(k) = config(k).opcoes(N);   % valor de preenchimento
                else
                    vars(k) = config(k).opcoes(idx);
                end

            elseif strcmp(modo_discreto, 'proporcional')
                % Heuristica legada (nao consta de [DF2011])
                n_estados = 2 ^ n;
                idx = floor(x_int * N / n_estados) + 1;
                idx = min(max(idx, 1), N);
                vars(k) = config(k).opcoes(idx);

            else
                error('decodificar:modoInvalido', ...
                    'modo_discreto deve ser ''datta'' ou ''proporcional''. Recebido: ''%s''.', ...
                    modo_discreto);
            end
        end
    end
end

end


% -------------------------------------------------------------------------
% >>> velocidade_binaria — Eq. (7)-(8), Tabelas 1-2 [DF2011]
% -------------------------------------------------------------------------

function [v_next, x_next] = velocidade_binaria(x_atual, tendencia, pm)
% RID_VELOCIDADE_BINARIA  Mapeia a velocidade real para {-1,0,+1} numa dimensao binaria.
%
% REFERENCIA: [DF2011] Secao 4.3, Eq. (7), Eq. (8), Tabelas 1 e 2.
%
% -------------------------------------------------------------------------
% Eq. (7) [DF2011]: transicoes de posicao binaria admissiveis
%
%   x^(i,t+1) = 0  se (x^(i,t), v^(i,t+1)) = (0,0) ou (1,-1)
%   x^(i,t+1) = 1  se (x^(i,t), v^(i,t+1)) = (0,1) ou (1,0)
%
%   Ou seja:   x = 0  =>  v pertence a { 0, +1}     (Tabela 1)
%              x = 1  =>  v pertence a { 0, -1}     (Tabela 2)
%
% [D3] Esta assimetria e o ponto central: a versao anterior deste solver
%      usava a regra SIMETRICA v em {-1,0,+1} independentemente do bit atual,
%      saturando a posicao depois. Isso corrompia a memoria de inercia
%      (w * v no passo seguinte) em varios casos das Tabelas 1 e 2.
%
% -------------------------------------------------------------------------
% Eq. (8) [DF2011]: mapeamento do valor real de velocidade para o discreto
%
%   v > 0                          =>  v = +1
%   v < 0                          =>  v = -1
%   v >= 0, ou sem relacao na Tab.1 =>  v = 0 ou +1
%   v <= 0, ou sem relacao na Tab.2 =>  v = 0 ou -1
%
% A ambiguidade ("0 ou 1" / "0 ou -1") e resolvida pela probabilidade de
% mutacao p_m (coluna 7 das Tabelas 1 e 2): com probabilidade p_m adota-se o
% valor ALTERNATIVO ao mapeado deterministicamente; caso contrario, o mapeado.
%
% [D4] A versao anterior sorteava uniformemente em {-1,0,+1} sem condicionar
%      ao bit atual, desperdicando parte dos sorteios na saturacao posterior.
%
% [DF2011] Sec. 4.3 sobre p_m: "Extensive empirical studies have shown that
% better results are obtained for p_m <= 15%."
%
% -------------------------------------------------------------------------
% ENTRADAS
%   x_atual   : bit atual da dimensao (0 ou 1)
%   tendencia : valor REAL da velocidade calculado pela Eq. (1)
%               v = w*v + c1*r1*(pbest - x) + c2*r2*(gbest - x)
%   pm        : probabilidade de mutacao (adocao do valor alternativo)
%
% SAIDAS
%   v_next : velocidade discreta resultante (-1, 0 ou +1)
%   x_next : nova posicao binaria, x_atual + v_next (garantidamente 0 ou 1)

validateattributes(x_atual, {'numeric'}, {'scalar'});
assert(x_atual == 0 || x_atual == 1, ...
    'velocidade_binaria:posicaoInvalida', ...
    'x_atual deve ser 0 ou 1. Recebido: %g.', x_atual);

if x_atual == 0
    % ---- Tabela 1 [DF2011]: valores admissiveis de v sao 0 ou +1 ----
    if tendencia > 0
        v_det = 1;
    elseif tendencia < 0
        v_det = 0;    % -1 e inadmissivel para x=0 (Eq. 7)
    else
        v_det = double(rand < 0.5);   % empate exato: 50% cada
    end
    v_alt = 1 - v_det;                % alternativa dentro de {0, +1}

else
    % ---- Tabela 2 [DF2011]: valores admissiveis de v sao 0 ou -1 ----
    if tendencia < 0
        v_det = -1;
    elseif tendencia > 0
        v_det = 0;    % +1 e inadmissivel para x=1 (Eq. 7)
    else
        v_det = -double(rand < 0.5);  % empate exato: 50% cada
    end
    v_alt = -1 - v_det;               % alternativa dentro de {0, -1}
end

% Randomizacao pela probabilidade de mutacao p_m
% (coluna 7 das Tabelas 1 e 2 [DF2011])
if rand < pm
    v_next = v_alt;
else
    v_next = v_det;
end

% Eq. (2)/(7) [DF2011]: por construcao das transicoes admissiveis, a soma
% resulta sempre em {0,1} — nao ha necessidade de saturacao.
x_next = x_atual + v_next;

end


% -------------------------------------------------------------------------
% >>> mutacao_polinomial — Eq. (9) [DF2011]
% -------------------------------------------------------------------------

function x_novo = mutacao_polinomial(x, x_l, x_u, eta)
% RID_MUTACAO_POLINOMIAL  Mutacao polinomial de um numero real.
%
% REFERENCIA: [DF2011] Eq. (9), Secao 5.
%   Originalmente de: Deb, K. "Multi-Objective Optimization using
%   Evolutionary Algorithms", John Wiley & Sons, 2001 (ref. [9] de [DF2011]).
%
% -------------------------------------------------------------------------
% Eq. (9) [DF2011]:
%
%   x  <-  x + (x^(u) - x^(l)) * d_q
%
%   d_q = [ 2r + (1-2r) * ( (x^(u) - x)/(x^(u) - x^(l)) )^(eta+1) ]^(1/(eta+1)) - 1
%         se r < 0.5
%
%   d_q = 1 - [ 2(1-r) + (2r-1) * ( (x - x^(l))/(x^(u) - x^(l)) )^(eta+1) ]^(1/(eta+1))
%         caso contrario
%
%   onde r e um numero aleatorio em ]0,1[ e eta > 0 o indice de distribuicao
%   polinomial.
%
% [DF2011] Sec. 5: "The distribution index eta in Eq. (9) is assigned a random
% value, in the range of [25,45], in different runs of a problem."
%
% [D5] A versao anterior deste solver usava perturbacao GAUSSIANA simples
%      (randn * 5% do dominio) na busca local, em vez da Eq. (9).
%
% -------------------------------------------------------------------------
% ENTRADAS
%   x   : valor real atual
%   x_l : limite inferior do dominio
%   x_u : limite superior do dominio
%   eta : indice de distribuicao polinomial (> 0)
%
% SAIDA
%   x_novo : valor mutado, garantidamente dentro de [x_l, x_u]
%
% PROPRIEDADE: quanto MAIOR eta, menor a perturbacao (distribuicao mais
% concentrada em torno de x). Essa propriedade e verificada em
% 06_testes/TestMutacaoPolinomial.m.

if x_u <= x_l
    % Dominio degenerado: nada a mutar
    x_novo = x;
    return;
end

assert(eta > 0, 'mutacao_polinomial:etaInvalido', ...
    'O indice de distribuicao eta deve ser positivo. Recebido: %g.', eta);

r     = rand;
delta = x_u - x_l;

if r < 0.5
    base = 2*r + (1 - 2*r) * ((x_u - x) / delta)^(eta + 1);
    d_q  = base^(1/(eta + 1)) - 1;
else
    base = 2*(1 - r) + (2*r - 1) * ((x - x_l) / delta)^(eta + 1);
    d_q  = 1 - base^(1/(eta + 1));
end

x_novo = x + delta * d_q;

% Salvaguarda numerica: a Eq. (9) e construida para manter x em [x_l, x_u],
% mas arredondamentos podem levar a excursoes de ordem 1e-16.
x_novo = max(x_l, min(x_u, x_novo));

end


% -------------------------------------------------------------------------
% >>> domina_deb — [DEB2000] pag. 316
% -------------------------------------------------------------------------

function melhor = domina_deb(custo_novo, viol_novo, custo_ref, viol_ref)
% RID_DOMINA_DEB  Regra de viabilidade (dominancia) de Deb, sem parametro de penalizacao.
%
% REFERENCIA:
%   [DEB2000] Deb, K. "An efficient constraint handling method for genetic
%             algorithms". Comput. Methods Appl. Mech. Engrg. 186 (2000)
%             311-338.  ->  pag. 316
%
% Citacao literal dos tres criterios ([DEB2000], pag. 316):
%   "1. Any feasible solution is preferred to any infeasible solution.
%    2. Among two feasible solutions, the one having better objective
%       function value is preferred.
%    3. Among two infeasible solutions, the one having smaller constraint
%       violation is preferred."
%
% Nao utiliza parametro de penalizacao (abordagem "penalty-parameter-less"),
% adotada tambem por [DF2011] Sec. 5: "the penalty-parameter-less constraint
% handling approach, proposed by Deb [8], is applied here for working with
% infeasible solutions."
%
% ENTRADAS
%   custo_novo, viol_novo : candidato
%   custo_ref,  viol_ref  : referencia (p-best, g-best, etc.)
%
% SAIDA
%   melhor : true se o candidato e PREFERIVEL a referencia
%
% CONVENCAO: violacao <= 0 significa solucao VIAVEL.

novo_viavel = (viol_novo <= 0);
ref_viavel  = (viol_ref  <= 0);

if novo_viavel && ~ref_viavel
    % Criterio 1: viavel domina inviavel
    melhor = true;

elseif ~novo_viavel && ref_viavel
    % Criterio 1: inviavel nunca domina viavel
    melhor = false;

elseif novo_viavel && ref_viavel
    % Criterio 2: ambos viaveis -> menor custo vence
    melhor = (custo_novo < custo_ref);

else
    % Criterio 3: ambos inviaveis -> menor violacao vence
    melhor = (viol_novo < viol_ref);
end

end


% -------------------------------------------------------------------------
% >>> amostrar_particula — [DF2011] Sec. 4.1
% -------------------------------------------------------------------------

function p = amostrar_particula(mapa, config, total_dim)
% Inicializacao aleatoria de uma particula. [DF2011] Sec. 4.1:
% dimensoes reais recebem valor real aleatorio no dominio; dimensoes
% binarias recebem 0 ou 1 com 50% de probabilidade.

p = zeros(1, total_dim);
for k = 1:length(mapa)
    info = mapa(k);
    orig = info.idx_original;
    if strcmp(info.tipo, 'R')
        vmin = config(orig).min;
        vmax = config(orig).max;
        p(info.inicio) = vmin + rand * (vmax - vmin);
    else
        p(info.inicio:info.fim) = double(rand(1, info.n_bits) < 0.5);
    end
end
end


% -------------------------------------------------------------------------
% >>> get_param — utilitario
% -------------------------------------------------------------------------

function val = get_param(s, field_name, default_val)
% Le um campo de hiperparametro do struct pso_params, com valor padrao.

if isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
else
    val = default_val;
end
end
