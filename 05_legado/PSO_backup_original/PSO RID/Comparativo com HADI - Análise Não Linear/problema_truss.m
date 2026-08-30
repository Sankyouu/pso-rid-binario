function problema_truss

% OTIMIZAÇÃO DA TRELIÇA DE 10 BARRAS — Hadi & Alvani (2003), Seção 6.1
% Case 2 (catálogo 2, variáveis discretas, análise linear)
% Resultado esperado (Tabela 1, Case 2, método proposto):
%   Áreas [mm²]: 19355 | 65 | 16129 | 7742 | 65 | 645 | 5161 | 12903 | 12903 | 65
%   Peso        : 2325.2 kg
%
% RESTRIÇÕES (Seção 6.1 do artigo):
%   Tensão admissível  : |sigma_i| ≤ 172.25 MPa  (tração e compressão)
%   Deslocamento máx.  : |u_j|    ≤  50.80 mm    (todos os DOFs livres)

clc; close all;

% --- CATÁLOGO DE SEÇÕES DISPONÍVEIS [mm²] ---
 catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, 11613, 12903, 16129, ...
            19355, 22581, 29032];

% --- LIMITES DAS RESTRIÇÕES ---
sigma_max = 172.25; % [MPa] — tensão admissível (tração e compressão)
d_max     =  50.80; % [mm]  — deslocamento admissível

lb = 65;   % Área mínima de cada barra [in²]
ub = 29032;  % Área máxima de cada barra [in²]
n_barras = 10

% --- CONFIG DO SOLVER: 10 variáveis discretas (uma por barra) ---

%for i = 1:n_barras
%    vars_solver(i).tipo = 'R';
%    vars_solver(i).min  = lb;
%    vars_solver(i).max  = ub;
%end

for i = 1:10
    vars_solver(i).tipo   = 'D';
    vars_solver(i).opcoes = catalogo;
end

%fprintf('============================================================\n');
%fprintf(' TRELIÇA 10 BARRAS — PSO-RID  (Hadi & Alvani, 2003)\n');
%fprintf(' Catálogo  : %d seções disponíveis\n', length(catalogo));
%fprintf(' sigma_max : %.2f MPa\n', sigma_max);
%fprintf(' d_max     : %.2f mm\n',  d_max);
%fprintf(' Ref. artigo (Case 2): 2325.2 kg\n');
%fprintf('============================================================\n\n');

% --- MULTI-START ---
n_runs             = 10;
melhor_peso_global = inf;
melhor_sol_global  = [];
todos_historicos   = {};

funcao_para_solver = @(x) custo_truss(x, sigma_max, d_max);

for run = 1:n_runs
    fprintf('--- Run %d/%d ---\n', run, n_runs);
    [sol, peso, hist] = pso_rid_generico(funcao_para_solver, vars_solver);
    todos_historicos{run} = hist;

    if peso < melhor_peso_global
        melhor_peso_global = peso;
        melhor_sol_global  = sol;
    end
    fprintf('  Run %d: peso = %.2f kg  |  areas = [', run, peso);
    fprintf(' %g', sol);
    fprintf(' ]\n\n');
end

% --- RESULTADOS FINAIS ---
[~, viol_final]                    = custo_truss(melhor_sol_global, sigma_max, d_max);
[weight_final, Sigma_final, u_final] = truss(melhor_sol_global);

% Deslocamentos máximos por DOF livre
n_dofs_livres = length(u_final); % 8 DOFs livres (12 - 4 apoios)

% Resultado de referência do artigo (Case 2, catálogo 2)
ref_areas = [19355, 65, 16129, 7742, 65, 645, 5161, 12903, 12903, 65];
ref_peso  = 2325.2;

fprintf('============================================================\n');
fprintf(' COMPARATIVO COM O ARTIGO\n');
fprintf('------------------------------------------------------------\n');
fprintf('         A1      A2      A3      A4      A5\n');
fprintf(' PSO : %6g  %6g  %6g  %6g  %6g\n', melhor_sol_global(1:5));
fprintf(' Ref : %6g  %6g  %6g  %6g  %6g\n', ref_areas(1:5));
fprintf('         A6      A7      A8      A9     A10\n');
fprintf(' PSO : %6g  %6g  %6g  %6g  %6g\n', melhor_sol_global(6:10));
fprintf(' Ref : %6g  %6g  %6g  %6g  %6g\n', ref_areas(6:10));
fprintf('------------------------------------------------------------\n');
fprintf(' Peso PSO    : %.2f kg\n', weight_final);
fprintf(' Peso artigo : %.2f kg\n', ref_peso);
fprintf(' Diferença   : %.2f kg  (%.2f%%)\n', ...
        weight_final - ref_peso, 100*(weight_final - ref_peso)/ref_peso);
fprintf(' Violação    : %.6e\n', viol_final);
fprintf('------------------------------------------------------------\n');
fprintf(' Barra  Área[mm²]   Tensão[MPa]   |sigma|/sigma_max\n');
fprintf('------------------------------------------------------------\n');
for i = 1:10
    fprintf('  %2d    %6g      %10.3f      %.3f\n', ...
            i, melhor_sol_global(i), Sigma_final(i), ...
            abs(Sigma_final(i))/sigma_max);
end
fprintf('------------------------------------------------------------\n');
fprintf(' DOF livre  Deslocamento[mm]   |u|/d_max\n');
fprintf('------------------------------------------------------------\n');
% DOFs livres: 1-8 correspondem aos nós 1..4 (x e y de cada)
% Mapeamento: DOF global → nó e direção
nos_livres = [1 1 2 2 3 3 4 4];   % nó de cada DOF livre
dir_livres = {'x','y','x','y','x','y','x','y'}; % direção
for j = 1:n_dofs_livres
    flag = '';
    if abs(u_final(j)) > d_max
        flag = '  *** VIOLA ***';
    end
    fprintf('  %2d (nó %d, %s)   %10.4f          %.4f%s\n', ...
            j, nos_livres(j), dir_livres{j}, u_final(j), ...
            abs(u_final(j))/d_max, flag);
end
fprintf(' Deslocamento máx. encontrado : %.4f mm  (limite: %.2f mm)\n', ...
        max(abs(u_final)), d_max);
fprintf('============================================================\n');

% --- PLOTS ---
figure('Name', 'PSO-RID — Treliça Hadi 2003', 'NumberTitle', 'off', ...
       'Position', [100 100 900 500]);

subplot(1,2,1);
cores = lines(n_runs);
hold on;
for run = 1:n_runs
    semilogy(todos_historicos{run}, 'Color', [cores(run,:), 0.6], 'LineWidth', 1.2);
end
yline(ref_peso, 'k--', 'LineWidth', 1.5);
xlabel('Iteração'); ylabel('Peso [kg] (escala log)');
title('Convergência por Run');
labels = arrayfun(@(x) sprintf('Run %d', x), 1:n_runs, 'UniformOutput', false);
labels{end+1} = 'Ref. artigo';
legend(labels, 'Location', 'northeast', 'FontSize', 6);
grid on; box on;

[~, idx_melhor] = min(cellfun(@(h) h(end), todos_historicos));
subplot(1,2,2);
semilogy(todos_historicos{idx_melhor}, 'b-', 'LineWidth', 2); hold on;
yline(ref_peso, 'r--', 'LineWidth', 1.5);
xlabel('Iteração'); ylabel('Peso [kg] (escala log)');
title(sprintf('Melhor Run (%d) — %.2f kg', idx_melhor, melhor_peso_global));
legend('PSO-RID', 'Ref. artigo', 'Location', 'northeast');
grid on; box on;

sgtitle('Treliça 10 Barras — Hadi & Alvani (2003)', ...
        'FontSize', 13, 'FontWeight', 'bold');
end


% FUNÇÃO DE CUSTO
function [J, violacao] = custo_truss(area, sigma_max, d_max)
% Retorna:
%   J        — Peso [kg] (função objetivo)
%   violacao — Soma ponderada das violações de tensão E deslocamento
%              (tratada via Método de Deb no solver)

[weight, Sigma, u_livre] = truss(area);

violacao = 0;

% --- Restrição de TENSÃO ---
% |sigma_i| ≤ sigma_max
for i = 1:length(Sigma)
    excesso = abs(Sigma(i)) - sigma_max;
    if excesso > 0
        violacao = violacao + excesso^2;
    end
end

% --- Restrição de DESLOCAMENTO ---
% |u_j| ≤ d_max
for j = 1:length(u_livre)
    excesso = abs(u_livre(j)) - d_max;
    if excesso > 0
        violacao = violacao + excesso^2;
    end
end

J = weight; % Peso em kg — o Método de Deb no solver trata restrições via violação
end
