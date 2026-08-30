function problema_truss

clc; close all;

% --- FATOR DE CONVERSÃO ---
in2_to_mm2 = 645.16; % 1 in² = 645.16 mm²

% --- CATÁLOGOS (convertidos de in² para mm²) ---
cat_A1 = [21.5, 22.5, 23.5, 24.5]             * in2_to_mm2;
cat_A2 = [0.1,  0.15, 0.2,  0.25]             * in2_to_mm2;
cat_A3 = [22.4, 25.4, 27.4, 29.4]             * in2_to_mm2;
cat_A4 = [14.1, 14.2, 14.3, 14.4, 14.5]       * in2_to_mm2;
cat_A6 = [0.5,  1.0,  1.5,  2.0,  2.5]        * in2_to_mm2;
cat_A7 = [11.0, 11.3, 11.7, 12.0, 12.3, 12.5] * in2_to_mm2;

% --- TENSÃO ADMISSÍVEL ---
sigma_max = 172.25; % [MPa]

% --- CONFIG DO SOLVER (uma struct por barra) ---

% A1 — discreta
config_vars(1).tipo   = 'D';
config_vars(1).opcoes = cat_A1;

% A2 — discreta
config_vars(2).tipo   = 'D';
config_vars(2).opcoes = cat_A2;

% A3 — discreta
config_vars(3).tipo   = 'D';
config_vars(3).opcoes = cat_A3;

% A4 — discreta
config_vars(4).tipo   = 'D';
config_vars(4).opcoes = cat_A4;

% A5 = A2 — discreta (mesmo catálogo)
config_vars(5).tipo   = 'D';
config_vars(5).opcoes = cat_A2;

% A6 — discreta
config_vars(6).tipo   = 'D';
config_vars(6).opcoes = cat_A6;

% A7 — discreta
config_vars(7).tipo   = 'D';
config_vars(7).opcoes = cat_A7;

% A8 = A7 — discreta (mesmo catálogo)
config_vars(8).tipo   = 'D';
config_vars(8).opcoes = cat_A7;

% A9 — real contínua entre 19 e 25 in² → convertido para mm²
config_vars(9).tipo = 'R';
config_vars(9).min  = 19.0 * in2_to_mm2;
config_vars(9).max  = 25.0 * in2_to_mm2;

% A10 = A2 — discreta (mesmo catálogo)
config_vars(10).tipo   = 'D';
config_vars(10).opcoes = cat_A2;

fprintf('============================================================\n');
fprintf(' TRELIÇA 10 BARRAS — Config mista por barra\n');
fprintf(' sigma_max : %.2f MPa\n', sigma_max);
fprintf('------------------------------------------------------------\n');
fprintf(' Barra  Tipo   Opções / Intervalo\n');
fprintf('------------------------------------------------------------\n');
nomes_tipo = {'D','D','D','D','D','D','D','D','R','D'};
for i = 1:10
    if strcmp(config_vars(i).tipo, 'R')
        fprintf('  A%-2d    R     [%.2f, %.2f] contínuo\n', ...
                i, config_vars(i).min, config_vars(i).max);
    else
        ops = config_vars(i).opcoes;
        str = sprintf('%.4g ', ops);
        fprintf('  A%-2d    D     {%s}\n', i, strtrim(str));
    end
end
fprintf('============================================================\n\n');

% --- MULTI-START ---
n_runs             = 10;
melhor_peso_global = inf;
melhor_sol_global  = [];
todos_historicos   = {};

funcao_custo = @(x) custo_truss(x, sigma_max);

for run = 1:n_runs
    fprintf('--- Run %d/%d ---\n', run, n_runs);
    [sol, peso, hist] = pso_rid_generico(funcao_custo, config_vars);
    todos_historicos{run} = hist;

    if peso < melhor_peso_global
        melhor_peso_global = peso;
        melhor_sol_global  = sol;
    end
    fprintf('  Run %d: peso = %.4f kg  |  areas = [', run, peso);
    fprintf(' %.4g', sol);
    fprintf(' ]\n\n');
end

% --- RESULTADOS FINAIS ---
[~, viol_final]             = custo_truss(melhor_sol_global, sigma_max);
[weight_final, Sigma_final] = truss(melhor_sol_global);

fprintf('============================================================\n');
fprintf(' MELHOR SOLUÇÃO ENCONTRADA\n');
fprintf('------------------------------------------------------------\n');
fprintf(' Peso total : %.4f kg\n', weight_final);
fprintf(' Violação   : %.6e\n',   viol_final);
fprintf('------------------------------------------------------------\n');
fprintf(' Barra  Tipo  Área[mm²]   Tensão[MPa]   |sigma|/sigma_max\n');
fprintf('------------------------------------------------------------\n');
for i = 1:10
    if strcmp(config_vars(i).tipo, 'R')
        tipo_str = 'R';
    else
        tipo_str = 'D';
    end
    fprintf('  A%-2d    %s    %7.4f     %10.4f       %.3f\n', ...
            i, tipo_str, melhor_sol_global(i), Sigma_final(i), ...
            abs(Sigma_final(i))/sigma_max);
end
fprintf('============================================================\n');

% --- PLOTS ---
figure('Name', 'PSO-RID — Treliça Config Mista', 'NumberTitle', 'off', ...
       'Position', [100 100 900 500]);

subplot(1,2,1);
cores = lines(n_runs);
hold on;
for run = 1:n_runs
    semilogy(todos_historicos{run}, 'Color', [cores(run,:), 0.6], 'LineWidth', 1.2);
end
xlabel('Iteração'); ylabel('Peso [kg] (escala log)');
title('Convergência por Run');
labels = arrayfun(@(r) sprintf('Run %d', r), 1:n_runs, 'UniformOutput', false);
legend(labels, 'Location', 'northeast', 'FontSize', 6);
grid on; box on;

[~, idx_melhor] = min(cellfun(@(h) h(end), todos_historicos));
subplot(1,2,2);
semilogy(todos_historicos{idx_melhor}, 'b-', 'LineWidth', 2);
xlabel('Iteração'); ylabel('Peso [kg] (escala log)');
title(sprintf('Melhor Run (%d) — %.4f kg', idx_melhor, melhor_peso_global));
legend('PSO-RID', 'Location', 'northeast');
grid on; box on;

sgtitle('Treliça 10 Barras — Config Mista por Barra', ...
        'FontSize', 13, 'FontWeight', 'bold');
end


% FUNÇÃO DE CUSTO
function [J, violacao] = custo_truss(area, sigma_max)

[weight, Sigma] = truss(area);

violacao = 0;
for i = 1:length(Sigma)
    excesso = abs(Sigma(i)) - sigma_max;
    if excesso > 0
        violacao = violacao + excesso^2;
    end
end

J = weight;
end