function run_awruch_discreto
%% OTIMIZAÇÃO DA TRELIÇA DE 10 BARRAS COM CATÁLOGO AWRUCH
% Otimização com múltiplos catálogos discretos independentes por barra

clc; close all;

root_dir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(root_dir);
setup_paths(false);

caso = catalogo_awruch();
sigma_max = caso.sigma_max;
d_max     = caso.d_max;

% Função de custo não linear
funcao_objetivo = @(areas) custo_awruch(areas, sigma_max, d_max);

% Configuração do PSO-RID
pso_params.n_particulas   = 100;
pso_params.max_iter       = 1000;
pso_params.tol_estagnacao = 200;
pso_params.verbose        = true;

n_runs = 5;
melhor_peso_global = inf;
melhor_sol_global  = [];
todos_historicos   = cell(n_runs, 1);

fprintf('============================================================\n');
fprintf(' INICIANDO OTIMIZAÇÃO PSO-RID: %s\n', caso.nome);
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

[weight_final, Sigma_final, u_final] = fem_truss_nonlinear(melhor_sol_global);
[~, viol_final] = custo_awruch(melhor_sol_global, sigma_max, d_max);

relatorio_comparativo('Treliça 10 Barras (Catálogo Awruch)', ...
                      melhor_sol_global, weight_final, ...
                      [], NaN, ...
                      Sigma_final, sigma_max, ...
                      u_final, d_max, viol_final);

fig = plot_convergencia(todos_historicos, NaN, 'PSO-RID — Treliça Catálogo Awruch');
salvar_figura(fig, 'convergencia_awruch_nao_linear');

end

function [J, violacao] = custo_awruch(areas, sigma_max, d_max)
    [weight, Sigma, u_livre] = fem_truss_nonlinear(areas);
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
