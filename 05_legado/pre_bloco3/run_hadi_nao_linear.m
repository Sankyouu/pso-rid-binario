function run_hadi_nao_linear(seed)
% OTIMIZAÇÃO DA TRELIÇA DE 10 BARRAS (ANÁLISE NÃO LINEAR) COM PSO-RID
% Benchmark: Hadi & Alvani (2003) - Case 2 (Variáveis Discretas)
% Análise Estrutural: Geometricamente Não Linear (Newton-Raphson KT = KE + KG)
%
% Uso:
%   run_hadi_nao_linear        %% usa semente padrão (42), reprodutível
%   run_hadi_nao_linear(123)   %% permite especificar outra semente

clc; close all;

if nargin < 1 || isempty(seed)
    seed = 42;
end
rng(seed);

% Garantir que as pastas do projeto estejam no path
root_dir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(root_dir);
setup_paths(false);

% --- 1. CARREGAR CONFIGURAÇÃO DO BENCHMARK ---
caso = catalogo_hadi;
sigma_max = caso.sigma_max;
d_max     = caso.d_max;

% --- 2. DEFINIÇÃO DA FUNÇÃO OBJETIVO COM RESTRIÇÕES ---
% Usa a função de custo centralizada em 01_src/utils/custo_trelica.m
% (elimina a duplicação que antes existia como função local neste arquivo)
funcao_objetivo = @(areas) custo_trelica(areas, sigma_max, d_max, 'nao_linear');

% --- 3. CONFIGURAÇÃO DO OTIMIZADOR PSO-RID ---
pso_params.n_particulas   = 100;
pso_params.max_iter       = 1000;
% Valores iniciais conforme Datta & Figueira (2011), Seção 5, pág. 3628.
% Com auto_adaptativo = true (padrão), cada parâmetro é sorteado a cada
% iteração no intervalo [0, valor_inicial].
pso_params.w              = 1.0;   %% inércia inicial      [DF2011] Sec. 5
pso_params.c1             = 1.0;   %% fator cognitivo      [DF2011] Sec. 5
pso_params.c2             = 2.0;   %% fator social         [DF2011] Sec. 5
pso_params.pm             = 0.15;  %% prob. de mutação     [DF2011] Sec. 4.3/5
pso_params.tol_estagnacao = 200;
pso_params.verbose        = true;
pso_params.print_interval = 100;

% --- 4. EXECUÇÃO MULTI-START ---
n_runs             = 5;
melhor_peso_global = inf;
melhor_sol_global  = [];
todos_historicos   = cell(n_runs, 1);

fprintf('============================================================\n');
fprintf(' INICIANDO OTIMIZAÇÃO PSO-RID: %s (NÃO LINEAR)\n', caso.nome);
fprintf(' Semente (rng): %d | População: %d partículas | Max Iter: %d | Runs: %d\n', ...
        seed, pso_params.n_particulas, pso_params.max_iter, n_runs);
fprintf('============================================================\n\n');

for r = 1:n_runs
    fprintf('--- Executando Run %d/%d ---\n', r, n_runs);
    [sol, peso, hist, details] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);
    todos_historicos{r} = hist;

    if peso < melhor_peso_global
        melhor_peso_global = peso;
        melhor_sol_global  = sol;
    end
    fprintf(' -> Fim do Run %d: Melhor Peso = %.2f kg (Iters: %d)\n\n', ...
            r, peso, details.iter_executadas);
end

% --- 5. AVALIAÇÃO DETALHADA DA MELHOR SOLUÇÃO ---
[weight_final, Sigma_final, u_final] = fem_truss_nonlinear(melhor_sol_global);
[~, viol_final] = custo_trelica(melhor_sol_global, sigma_max, d_max, 'nao_linear');

% --- 6. EXIBIÇÃO DO RELATÓRIO COMPARATIVO ---
relatorio_comparativo('Treliça 10 Barras Hadi (Não Linear)', ...
                      melhor_sol_global, weight_final, ...
                      caso.ref_areas, caso.ref_peso, ...
                      Sigma_final, sigma_max, ...
                      u_final, d_max, viol_final);

% --- 7. PLOT E SALVAMENTO DOS RESULTADOS ---
fig = plot_convergencia(todos_historicos, caso.ref_peso, ...
                        'PSO-RID — Treliça 10 Barras (Análise Não Linear)');
salvar_figura(fig, 'convergencia_hadi_nao_linear');

end
