function [best_sol, best_cost, cost_history, details] = pso_rid(funcao_custo, config_vars, pso_params)
% PSO_RID  Particle Swarm Optimization Real-Inteiro-Discreto (codificacao binaria).
%
% BLOCO 1 do projeto CPIO III. Este arquivo contem apenas o LOOP do PSO;
% cada peca do algoritmo esta isolada em um auxiliar proprio (testavel
% individualmente por matlab.unittest):
%
%   rid_mapear_dimensoes.m     Eq. (5)      dimensoes da particula
%   rid_decodificar.m          Eq. (6)      binario -> real/inteiro/discreto
%   rid_velocidade_binaria.m   Eq. (7)-(8)  velocidade em dimensao binaria
%   rid_mutacao_polinomial.m   Eq. (9)      mutacao polinomial (variavel real)
%   rid_domina_deb.m           pag. 316     regra de viabilidade de Deb
%
% =========================================================================
% REFERENCIAS TEORICAS
% =========================================================================
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
% Mapeamento equacao -> local:
%   Eq. (1) [DF2011] velocidade real ............... secao C
%   Eq. (2) [DF2011] posicao real .................. secao C
%   Eq. (3) [DF2011] atualizacao do p-best ......... secao A
%   Eq. (4) [DF2011] atualizacao do g-best ......... secao A
%   Eq. (5) [DF2011] dimensoes da particula ........ rid_mapear_dimensoes
%   Eq. (6) [DF2011] decodificacao binaria ......... rid_decodificar
%   Eq. (7)-(8), Tab. 1-2 [DF2011] veloc. binaria .. rid_velocidade_binaria
%   Eq. (9) [DF2011] mutacao polinomial ............ rid_mutacao_polinomial
%   Sec. 5  [DF2011] busca local no g-best ......... secao B
%   pag. 316 [DEB2000] dominancia (3 criterios) .... rid_domina_deb
%
% =========================================================================
% ENTRADAS
% =========================================================================
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
%   cost_history : historico de convergencia do g-best (vetor coluna)
%   details      : struct diagnostico
%
% =========================================================================
% NOTA DE FIDELIDADE (revisao de 2026-08-28)
% =========================================================================
% A versao anterior (05_legado/pso_rid_pre_datta_2026_08_28.m) divergia de
% [DF2011] em cinco pontos, corrigidos aqui e marcados como [D1]..[D5] nos
% auxiliares correspondentes. Ver 00_docs/notas_e_relatorios/.

% =========================================================================
% VALIDACAO DE ENTRADA
% =========================================================================
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
                       'decodificacao_discreta', 'verbose', 'print_interval'};
campos_desconhecidos = setdiff(fieldnames(pso_params), campos_reconhecidos);
if ~isempty(campos_desconhecidos)
    warning('pso_rid:parametroDesconhecido', ...
        ['Campo(s) de pso_params nao reconhecido(s) e IGNORADO(S): %s. ' ...
         'Campos validos: %s.'], ...
        strjoin(campos_desconhecidos, ', '), strjoin(campos_reconhecidos, ', '));
end

% =========================================================================
% HIPERPARAMETROS
% -------------------------------------------------------------------------
% [DF2011] Sec. 5, pag. 3628: "In each run, the mutation probability p_m for
% binary variables is assigned the initial value of 15%. For handling real
% variables, the assigned initial values of the inertia constant w, cognitive
% behavioral factor c1, and social behavioral factor c2 are 1.0, 1.0 and 2.0,
% respectively."
% =========================================================================
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

% Estrategia de decodificacao discreta — ver [D1] em rid_decodificar.m
%   'datta'        (padrao) fiel a [DF2011] Sec. 4.1
%   'proporcional' heuristica legada, sem desperdicio de codigos
decod_discreta  = get_param(pso_params, 'decodificacao_discreta', 'datta');
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

% =========================================================================
% MAPEAMENTO DA PARTICULA — Eq. (5) [DF2011]
% =========================================================================
[mapa_dimensoes, total_dim] = rid_mapear_dimensoes(config_vars);

% Lookup O(1) de tipo / variavel de origem por dimensao
tipo_por_dim = cell(total_dim, 1);
orig_por_dim = zeros(total_dim, 1);
for k = 1:length(mapa_dimensoes)
    for d = mapa_dimensoes(k).inicio : mapa_dimensoes(k).fim
        tipo_por_dim{d} = mapa_dimensoes(k).tipo;
        orig_por_dim(d) = mapa_dimensoes(k).idx_original;
    end
end

% =========================================================================
% INICIALIZACAO DO ENXAME
% -------------------------------------------------------------------------
% [DF2011] Sec. 4.1: "A real dimension is initialized by a random real value
% in the given range of the real variable which is represented by that
% dimension, while all the binary dimensions are initialized randomly by 0 or
% 1 with 50% probability. On the other hand, all the velocity components are
% assigned the initial value of 0."
% =========================================================================
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
iter_sem_melhora = 0;
melhor_custo_ant = inf;
iter_final       = max_iter;
n_avaliacoes     = 0;
n_overflow       = 0;

% =========================================================================
% LOOP PRINCIPAL
% =========================================================================
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

    % =====================================================================
    % A. AVALIACAO E ATUALIZACAO DE p-BEST / g-BEST
    %    Eq. (3) e (4) [DF2011], dominancia de [DEB2000] pag. 316
    % =====================================================================
    for i = 1:n_particulas
        [vars_i, viol_estrutural] = rid_decodificar(pos(i,:), mapa_dimensoes, ...
                                                    config_vars, decod_discreta);

        if viol_estrutural > 0
            % Valor fora do dominio (ex.: indice inexistente no catalogo):
            % a solucao nao pode ser avaliada fisicamente. Tratada como
            % inviavel, sem calcular a funcao objetivo — coerente com
            % [DEB2000] pag. 316.
            custo_atual = inf;
            viol_atual  = viol_estrutural;
            n_overflow  = n_overflow + 1;
        else
            [custo_atual, viol_atual] = funcao_custo(vars_i);
            n_avaliacoes = n_avaliacoes + 1;
        end

        % Eq. (3) [DF2011]: atualizacao do p-best
        if rid_domina_deb(custo_atual, viol_atual, pbest_custo(i), pbest_viol(i))
            pbest_custo(i) = custo_atual;
            pbest_viol(i)  = viol_atual;
            pbest_pos(i,:) = pos(i,:);
        end

        % Eq. (4) [DF2011]: atualizacao do g-best
        if rid_domina_deb(pbest_custo(i), pbest_viol(i), gbest_custo, gbest_viol)
            gbest_custo = pbest_custo(i);
            gbest_viol  = pbest_viol(i);
            gbest_pos   = pbest_pos(i,:);
        end
    end

    % =====================================================================
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
    % =====================================================================
    if ~isinf(gbest_custo) || ~isinf(gbest_viol)
        gbest_temp = gbest_pos;
        p_flip     = rand * p_busca_local;    % probabilidade em ]0, 0.05[

        for d = 1:total_dim
            orig_d = orig_por_dim(d);
            if strcmp(tipo_por_dim{d}, 'R')
                gbest_temp(d) = rid_mutacao_polinomial(gbest_temp(d), ...
                                    config_vars(orig_d).min, ...
                                    config_vars(orig_d).max, eta_run);
            else
                if rand < p_flip
                    gbest_temp(d) = 1 - gbest_temp(d);
                end
            end
        end

        [vars_mut, viol_mut_estr] = rid_decodificar(gbest_temp, mapa_dimensoes, ...
                                                    config_vars, decod_discreta);
        if viol_mut_estr > 0
            c_mut = inf;
            v_mut = viol_mut_estr;
            n_overflow = n_overflow + 1;
        else
            [c_mut, v_mut] = funcao_custo(vars_mut);
            n_avaliacoes = n_avaliacoes + 1;
        end

        if rid_domina_deb(c_mut, v_mut, gbest_custo, gbest_viol)
            gbest_pos   = gbest_temp;
            gbest_custo = c_mut;
            gbest_viol  = v_mut;
        end
    end

    cost_history(iter) = gbest_custo;

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
        iter_final   = iter;
        break;
    end

    % =====================================================================
    % C. ATUALIZACAO DE VELOCIDADE E POSICAO
    % =====================================================================
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

                % Limite de velocidade (metade do dominio): estabilizacao
                % numerica, nao especificada em [DF2011].
                v_lim    = (vmax_d - vmin_d) / 2;
                vel(i,d) = max(-v_lim, min(v_lim, tendencia));

                % Eq. (2) [DF2011]: x^(t+1) = x^(t) + v^(t+1)
                pos(i,d) = pos(i,d) + vel(i,d);
                pos(i,d) = max(vmin_d, min(vmax_d, pos(i,d)));

            else
                % ---------- VARIAVEL BINARIA ----------
                % [D3][D4] Eq. (7)-(8) e Tabelas 1-2 [DF2011]
                [v_next, x_next] = rid_velocidade_binaria(pos(i,d), tendencia, pm_atual);
                vel(i,d) = v_next;
                pos(i,d) = x_next;
            end
        end
    end
end

% =========================================================================
% SAIDAS
% =========================================================================
best_cost = gbest_custo;
[best_sol, viol_final_estrutural] = rid_decodificar(gbest_pos, mapa_dimensoes, ...
                                                    config_vars, decod_discreta);

details.iter_executadas = iter_final;
details.gbest_viol      = gbest_viol;
details.gbest_raw_pos   = gbest_pos;
details.mapa_dimensoes  = mapa_dimensoes;
details.total_dim       = total_dim;
details.n_avaliacoes    = n_avaliacoes;
details.n_overflow      = n_overflow;
details.eta_run         = eta_run;
details.viol_estrutural = viol_final_estrutural;
details.decodificacao   = decod_discreta;

end


% #########################################################################
% FUNCAO LOCAL (usada apenas por este arquivo)
% #########################################################################

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


function val = get_param(s, field_name, default_val)
if isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
else
    val = default_val;
end
end
