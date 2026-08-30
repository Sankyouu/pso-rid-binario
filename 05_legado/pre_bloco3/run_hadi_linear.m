function run_hadi_linear
%% OTIMIZAÇÃO DA TRELIÇA DE 10 BARRAS (ANÁLISE LINEAR) COM PSO-RID
% Benchmark: Hadi & Alvani (2003) - Case 2 (Variáveis Discretas)
% Análise Estrutural: Elástica Linear

clc; close all;

% Adiciona as pastas necessárias ao path
root_dir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(root_dir);
setup_paths(false);

% Carregar parâmetros
caso = catalogo_hadi();
sigma_max = caso.sigma_max;
d_max     = caso.d_max;

% Função de custo linear
funcao_objetivo = @(areas) custo_trelica_linear(areas, sigma_max, d_max);

% Configuração do PSO-RID
pso_params.n_particulas   = 100;
pso_params.max_iter       = 1000;
pso_params.tol_estagnacao = 200;
pso_params.verbose        = true;

% Execução
n_runs = 5;
melhor_peso_global = inf;
melhor_sol_global  = [];
todos_historicos   = cell(n_runs, 1);

fprintf('============================================================\n');
fprintf(' INICIANDO OTIMIZAÇÃO PSO-RID: %s (LINEAR)\n', caso.nome);
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

% Avaliação Final
[weight_final, Sigma_final, u_final] = fem_truss_linear(melhor_sol_global);
[~, viol_final] = custo_trelica_linear(melhor_sol_global, sigma_max, d_max);

% Relatório
relatorio_comparativo('Treliça 10 Barras Hadi (Linear)', ...
                      melhor_sol_global, weight_final, ...
                      caso.ref_areas, caso.ref_peso, ...
                      Sigma_final, sigma_max, ...
                      u_final, d_max, viol_final);

% Gráfico
fig = plot_convergencia(todos_historicos, caso.ref_peso, ...
                        'PSO-RID — Treliça 10 Barras (Análise Linear)');
salvar_figura(fig, 'convergencia_hadi_linear');

end

function [J, violacao] = custo_trelica_linear(areas, sigma_max, d_max)
    [weight, Sigma, u_livre] = fem_truss_linear(areas);
    violacao = 0;
    for i = 1:length(Sigma)
        exc = abs(Sigma(i)) - sigma_max;
        if exc > 0, violacao = violacao + exc^2; end
    end
    for j = 1:length(u_livre)
        exc = abs(u_livre(j)) - d_max;
        if exc > 0, violacao = violacao + exc^2; end
    end
    J = weight;
end
