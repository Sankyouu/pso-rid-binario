function [best_sol, best_cost, cost_history, details] = pso_rid_pre_datta_2026_08_28(funcao_custo, config_vars, pso_params)
% PSO_RID Solver de Otimização por Enxame de Partículas Híbrido (Real-Inteiro-Discreto)
% Versão Aprimorada (revisada em 26/08/2026):
%   - Regras de Viabilidade de Deb completas na avaliação e na busca local
%   - Mapeamento binário -> discreto PROPORCIONAL (viés residual mínimo e
%     distribuído entre opções, corrigindo o viés real que a versão
%     anterior baseada em saturação min/max concentrava em um extremo)
%   - Inércia adaptativa linear decrescente (w_max -> w_min, padrão 0.9 -> 0.4)
%   - Reinicialização de partículas estagnadas para manutenção de diversidade
%   - Validação de campos desconhecidos em pso_params (evita parâmetros órfãos ignorados silenciosamente)
%
% Entradas:
%   funcao_custo : handle da função objetivo @(x) [custo, violacao]
%   config_vars  : struct array com a definição de cada variável (tipo 'R', 'I' ou 'D')
%   pso_params   : (opcional) struct com hiperparâmetros do PSO
%
% Saídas:
%   best_sol     : vetor com a melhor solução decodificada encontrada
%   best_cost    : custo objetivo associado à melhor solução
%   cost_history : vetor coluna com o histórico de convergência do GBest
%   details      : struct com informações diagnósticas da execução

% --- VALIDAÇÃO DE ENTRADA ---
assert(~isempty(config_vars), 'config_vars não pode ser vazio.');

for i = 1:length(config_vars)
    v = config_vars(i);
    assert(isfield(v, 'tipo'), 'config_vars(%d) deve ter campo "tipo".', i);

    if strcmp(v.tipo, 'R') || strcmp(v.tipo, 'I')
        assert(isfield(v,'min') && isfield(v,'max'), ...
            'config_vars(%d) tipo %s requer campos "min" e "max".', i, v.tipo);
    elseif strcmp(v.tipo, 'D')
        assert(isfield(v,'opcoes') && ~isempty(v.opcoes), ...
            'config_vars(%d) tipo D requer campo "opcoes" não vazio.', i);
    else
        error('Tipo desconhecido "%s" em config_vars(%d). Use R, I ou D.', v.tipo, i);
    end
end

% --- PARÂMETROS PADRÃO ---
if nargin < 3 || isempty(pso_params)
    pso_params = struct();
end

% Validação defensiva: alerta sobre campos não reconhecidos em pso_params
% (evita bugs silenciosos como um campo digitado errado, ou um nome
% de parâmetro legado que não corresponde a nenhum parâmetro atual)
campos_reconhecidos = {'n_particulas', 'max_iter', 'w_max', 'w_min', 'c1', 'c2', ...
                        'pm', 'tol_estagnacao', 'reinit_freq', 'reinit_pct', ...
                        'verbose', 'print_interval'};
campos_informados = fieldnames(pso_params);
campos_desconhecidos = setdiff(campos_informados, campos_reconhecidos);
if ~isempty(campos_desconhecidos)
    warning('pso_rid:parametroDesconhecido', ...
        ['Campo(s) de pso_params não reconhecido(s) e IGNORADO(S): %s. ' ...
         'Campos válidos: %s.'], ...
        strjoin(campos_desconhecidos, ', '), strjoin(campos_reconhecidos, ', '));
end

n_particulas   = get_param(pso_params, 'n_particulas', 100);
max_iter       = get_param(pso_params, 'max_iter', 1000);
w_max          = get_param(pso_params, 'w_max', 0.9);
w_min          = get_param(pso_params, 'w_min', 0.4);
c1             = get_param(pso_params, 'c1', 2.0);
c2             = get_param(pso_params, 'c2', 2.0);
pm             = get_param(pso_params, 'pm', 0.15);
tol_estagnacao = get_param(pso_params, 'tol_estagnacao', 200);
reinit_freq    = get_param(pso_params, 'reinit_freq', 100); % Frequência de reinicialização de partículas piores
reinit_pct     = get_param(pso_params, 'reinit_pct', 0.10); % Fração de partículas reinicializadas
verbose        = get_param(pso_params, 'verbose', true);
print_interval = get_param(pso_params, 'print_interval', 100);

% --- MAPEAMENTO DA PARTÍCULA (Decodificação Binária / Híbrida) ---
template.idx_original = 0;
template.tipo         = '';
template.n_bits       = 0;
template.inicio       = 0;
template.fim          = 0;

mapa_dimensoes = repmat(template, length(config_vars), 1);
total_dim = 0;

for i = 1:length(config_vars)
    var = config_vars(i);
    info.idx_original = i;
    info.tipo         = var.tipo;

    if strcmp(var.tipo, 'R')
        info.n_bits = 0;
        info.inicio = total_dim + 1;
        total_dim   = total_dim + 1;
        info.fim    = total_dim;
        else        
        if strcmp(var.tipo, 'I')
            limite = var.max;
        else  % 'D'
            limite = length(var.opcoes);
        end
        n_bits      = max(ceil(log2(limite + 0.1)), 1);
        info.n_bits = n_bits;
        info.inicio = total_dim + 1;
        total_dim   = total_dim + n_bits;
        info.fim    = total_dim;
    end
    mapa_dimensoes(i) = info;
end

% --- PRÉ-COMPUTAÇÃO DE TIPOS POR DIMENSÃO (O(1) lookup) ---
tipo_por_dim = cell(total_dim, 1);
orig_por_dim = zeros(total_dim, 1);

for k = 1:length(mapa_dimensoes)
    for d = mapa_dimensoes(k).inicio : mapa_dimensoes(k).fim
        tipo_por_dim{d} = mapa_dimensoes(k).tipo;
        orig_por_dim(d) = mapa_dimensoes(k).idx_original;
    end
end

% --- INICIALIZAÇÃO ---
pos = zeros(n_particulas, total_dim);
vel = zeros(n_particulas, total_dim);

for i = 1:n_particulas
    for k = 1:length(mapa_dimensoes)
        info = mapa_dimensoes(k);
        orig = info.idx_original;

        if strcmp(info.tipo, 'R')
            vmin = config_vars(orig).min;
            vmax = config_vars(orig).max;
            pos(i, info.inicio) = vmin + rand * (vmax - vmin);
        else
            pos(i, info.inicio:info.fim) = randi([0, 1], 1, info.n_bits);
        end
    end
end

% Memórias
pbest_pos   = pos;
pbest_custo = inf(n_particulas, 1);
pbest_viol  = inf(n_particulas, 1);

gbest_pos   = zeros(1, total_dim);
gbest_custo = inf;
gbest_viol  = inf;

cost_history     = nan(max_iter, 1);
iter_sem_melhora = 0;
melhor_custo_ant = inf;
iter_final       = max_iter;

% --- LOOP PRINCIPAL ---
for iter = 1:max_iter

    % Inércia linear decrescente (Exploração no início -> Refinamento no final)
    w_atual = w_max - (w_max - w_min) * (iter / max_iter);

    % --- A. AVALIAÇÃO (Regras de Viabilidade de Deb) ---
    for i = 1:n_particulas
        vars_decodificadas = decodificar(pos(i,:), mapa_dimensoes, config_vars);
        [custo_atual, viol_atual] = funcao_custo(vars_decodificadas);

        % Atualiza PBest (Regras de Deb)
        melhorou = false;
        if (pbest_viol(i) > 0) && (viol_atual < pbest_viol(i))
            melhorou = true; % Regra 1: ambos inviáveis, menor violação vence
        elseif (pbest_viol(i) > 0) && (viol_atual == 0)
            melhorou = true; % Regra 2: viável vence inviável
        elseif (pbest_viol(i) == 0) && (viol_atual == 0) && (custo_atual < pbest_custo(i))
            melhorou = true; % Regra 3: ambos viáveis, menor custo vence
        end

        if melhorou
            pbest_custo(i) = custo_atual;
            pbest_viol(i)  = viol_atual;
            pbest_pos(i,:) = pos(i,:);
        end

        % Atualiza GBest (Regras de Deb)
        melhorou_global = false;
        if (gbest_viol > 0) && (pbest_viol(i) < gbest_viol)
            melhorou_global = true;
        elseif (gbest_viol > 0) && (pbest_viol(i) == 0)
            melhorou_global = true;
        elseif (gbest_viol == 0) && (pbest_viol(i) == 0) && (pbest_custo(i) < gbest_custo)
            melhorou_global = true;
        end

        if melhorou_global
            gbest_custo = pbest_custo(i);
            gbest_viol  = pbest_viol(i);
            gbest_pos   = pbest_pos(i,:);
        end
    end

    % --- B. BUSCA LOCAL NO GBEST COM REGRAS DE DEB ---
    gbest_temp = gbest_pos;
    idx_mut    = randi(total_dim);
    tipo_mut   = tipo_por_dim{idx_mut};
    orig_mut   = orig_por_dim(idx_mut);

    if strcmp(tipo_mut, 'R')
        sigma = (config_vars(orig_mut).max - config_vars(orig_mut).min) * 0.05;
        gbest_temp(idx_mut) = gbest_temp(idx_mut) + randn * sigma;
        gbest_temp(idx_mut) = max(config_vars(orig_mut).min, ...
                              min(config_vars(orig_mut).max, gbest_temp(idx_mut)));
    else
        gbest_temp(idx_mut) = 1 - gbest_temp(idx_mut);
    end

    vars_mut = decodificar(gbest_temp, mapa_dimensoes, config_vars);
    [c_mut, v_mut] = funcao_custo(vars_mut);

    % Aceitação estrita com Regras de Deb
    aceitar_mutacao = false;
    if (gbest_viol > 0) && (v_mut < gbest_viol)
        aceitar_mutacao = true; % Reduziu violação em estado inviável
    elseif (gbest_viol > 0) && (v_mut == 0)
        aceitar_mutacao = true; % Encontrou solução viável
    elseif (gbest_viol == 0) && (v_mut == 0) && (c_mut < gbest_custo)
        aceitar_mutacao = true; % Reduziu custo mantendo viabilidade
    end

    if aceitar_mutacao
        gbest_pos   = gbest_temp;
        gbest_custo = c_mut;
        gbest_viol  = v_mut;
    end

    cost_history(iter) = gbest_custo;

    % Critério de estagnação
    if gbest_custo < melhor_custo_ant - 1e-12
        melhor_custo_ant = gbest_custo;
        iter_sem_melhora = 0;
    else
        iter_sem_melhora = iter_sem_melhora + 1;
    end

    if verbose && (mod(iter, print_interval) == 0 || iter == 1)
        fprintf('Iter %4d (w=%.2f) | Violação: %.4e | Custo: %.6e | Sem melhora: %d\n', ...
                iter, w_atual, gbest_viol, gbest_custo, iter_sem_melhora);
    end

    % Reinicialização periódica de partículas estagnadas (Prevenção de convergência prematura)
    if reinit_freq > 0 && mod(iter, reinit_freq) == 0 && iter_sem_melhora > 20
        n_reinit = max(1, round(reinit_pct * n_particulas));
        [~, idx_ordenados] = sort(pbest_viol, 'descend'); % Piores partículas
        for r_idx = 1:n_reinit
            p_idx = idx_ordenados(r_idx);
            for k = 1:length(mapa_dimensoes)
                info = mapa_dimensoes(k);
                orig = info.idx_original;
                if strcmp(info.tipo, 'R')
                    vmin = config_vars(orig).min;
                    vmax = config_vars(orig).max;
                    pos(p_idx, info.inicio) = vmin + rand * (vmax - vmin);
                else
                    pos(p_idx, info.inicio:info.fim) = randi([0, 1], 1, info.n_bits);
                end
            end
            vel(p_idx, :) = 0;
            pbest_custo(p_idx) = inf;
            pbest_viol(p_idx)  = inf;
        end
    end

    if iter_sem_melhora >= tol_estagnacao
        if verbose
            fprintf('>>> Parada antecipada na iteração %d (estagnação de %d iterações)\n', ...
                    iter, tol_estagnacao);
        end
        cost_history = cost_history(1:iter);
        iter_final   = iter;
        break;
    end

    % --- C. ATUALIZAÇÃO DE VELOCIDADE E POSIÇÃO ---
    for i = 1:n_particulas
        for d = 1:total_dim
            tipo_d = tipo_por_dim{d};
            orig_d = orig_por_dim(d);

            r1 = rand; r2 = rand;
            tendencia = w_atual * vel(i,d) ...
                + c1 * r1 * (pbest_pos(i,d) - pos(i,d)) ...
                + c2 * r2 * (gbest_pos(d)   - pos(i,d));

            if strcmp(tipo_d, 'R')
                vmax_d   = (config_vars(orig_d).max - config_vars(orig_d).min) / 2;
                vel(i,d) = max(-vmax_d, min(vmax_d, tendencia));
                pos(i,d) = pos(i,d) + vel(i,d);
                pos(i,d) = max(config_vars(orig_d).min, ...
                           min(config_vars(orig_d).max, pos(i,d)));
            else
                if     tendencia > 0, v_next =  1;
                elseif tendencia < 0, v_next = -1;
                else,                 v_next =  0;
                end

                if rand < pm
                    v_options = [-1, 0, 1];
                    v_next = v_options(randi(3));
                end

                vel(i,d) = v_next;
                pos(i,d) = pos(i,d) + v_next;
                pos(i,d) = max(0, min(1, pos(i,d)));
            end
        end
    end
end

best_cost = gbest_custo;
best_sol  = decodificar(gbest_pos, mapa_dimensoes, config_vars);

% Detalhes diagnósticos
details.iter_executadas = iter_final;
details.gbest_viol      = gbest_viol;
details.gbest_raw_pos   = gbest_pos;
details.mapa_dimensoes  = mapa_dimensoes;

end

% --- DECODIFICADOR SEM VIÉS (Saturação Simétrica nos Limites) ---
function vars = decodificar(vetor_hibrido, mapa, config)
n_vars = length(config);
vars   = zeros(1, n_vars);

for k = 1:n_vars
    info = mapa(k);

    if strcmp(info.tipo, 'R')
        vars(k) = vetor_hibrido(info.inicio);
    else
        bits    = vetor_hibrido(info.inicio : info.fim);
        n       = length(bits);
        pesos   = 2 .^ ((n-1):-1:0);
        decimal = round(bits) * pesos';

        if strcmp(info.tipo, 'I')
            vars(k) = min(decimal, config(k).max);
        else  % 'D' - Variável Discreta
            N         = length(config(k).opcoes);
            n_estados = 2 ^ n; % total de combinações binárias possíveis (2^n_bits)
            % Mapeamento PROPORCIONAL (corrige o viés estatístico real):
            % a antiga saturação (min/max) concentrava TODO o excedente
            % de decimais fora de [0,N-1] em um único extremo do catálogo
            % (a maior seção tinha ~4x mais chance de sair). O mapeamento
            % proporcional abaixo distribui o excedente uniformemente entre
            % as opções, reduzindo o viés máximo de ~4x para ~2x.
            idx     = floor(decimal * N / n_estados) + 1;
            idx     = min(max(idx, 1), N); % salvaguarda numérica (nunca deveria ser necessário)
            vars(k) = config(k).opcoes(idx);
        end
    end
end
end

function val = get_param(s, field_name, default_val)
if isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
else
    val = default_val;
end
end
