function run_estudo_particulas
%% ESTUDO COMPARATIVO DO TAMANHO DO ENXAME (50 vs 100 PARTÍCULAS)
% Avalia o impacto do número de partículas na convergência e no peso final

clc; close all;

root_dir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(root_dir);
setup_paths(false);

caso = catalogo_hadi();
funcao_objetivo = @(areas) custo_trelica_nao_linear(areas, caso.sigma_max, caso.d_max);

tamanhos_pop = [50, 100];
n_runs = 5;
resultados = struct();

for p = 1:length(tamanhos_pop)
    n_part = tamanhos_pop(p);
    fprintf('\n============================================================\n');
    fprintf(' TESTANDO ENXAME COM %d PARTÍCULAS (%d RUNS)\n', n_part, n_runs);
    fprintf('============================================================\n');

    pso_params.n_particulas   = n_part;
    pso_params.max_iter       = 1000;
    pso_params.tol_estagnacao = 200;
    pso_params.verbose        = false;

    historicos = cell(n_runs, 1);
    pesos      = zeros(n_runs, 1);

    for r = 1:n_runs
        fprintf('  Executando Run %d/%d (Pop: %d)...', r, n_runs, n_part);
        [sol, peso, hist, details] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);
        historicos{r} = hist;
        pesos(r)      = peso;
        fprintf(' Peso = %.2f kg (Iters: %d)\n', peso, details.iter_executadas);
    end

    resultados(p).n_part     = n_part;
    resultados(p).historicos = historicos;
    resultados(p).pesos      = pesos;
    resultados(p).melhor     = min(pesos);
    resultados(p).media      = mean(pesos);
    resultados(p).desvio     = std(pesos);
end

% Comparação Estatística
fprintf('\n============================================================\n');
fprintf(' RESUMO COMPARATIVO: 50 vs 100 PARTÍCULAS\n');
fprintf('============================================================\n');
fprintf(' População | Melhor [kg] | Média [kg] | Desvio Padrão [kg]\n');
fprintf(' -----------------------------------------------------------\n');
for p = 1:length(resultados)
    fprintf('   %3d     |   %8.2f  |  %8.2f  |     %8.2f\n', ...
            resultados(p).n_part, resultados(p).melhor, ...
            resultados(p).media, resultados(p).desvio);
end
fprintf('============================================================\n\n');

% Gráfico Comparativo
fig = figure('Name', 'Estudo de Particulas', 'NumberTitle', 'off', 'Position', [150 150 900 450]);
cores = {'r', 'b'};

for p = 1:length(resultados)
    subplot(1, 2, p);
    hold on;
    for r = 1:n_runs
        semilogy(resultados(p).historicos{r}, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0);
    end
    [~, idx_melhor] = min(resultados(p).pesos);
    semilogy(resultados(p).historicos{idx_melhor}, cores{p}, 'LineWidth', 2.0);
    yline(caso.ref_peso, 'k--', 'LineWidth', 1.5);
    title(sprintf('%d Partículas (Melhor: %.2f kg)', resultados(p).n_part, resultados(p).melhor), 'FontWeight', 'bold');
    xlabel('Iteração'); ylabel('Peso [kg] (log)');
    legend('Runs', '', '', '', '', sprintf('Melhor (%d part)', resultados(p).n_part), 'Ref. Hadi', 'Location', 'northeast');
    grid on; box on;
end
sgtitle('Sensibilidade do Tamanho do Enxame no PSO-RID', 'FontWeight', 'bold');

salvar_figura(fig, 'estudo_sensibilidade_particulas');

end

function [J, violacao] = custo_trelica_nao_linear(areas, sigma_max, d_max)
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
