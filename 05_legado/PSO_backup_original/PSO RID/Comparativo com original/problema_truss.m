function problema_truss

clc; close all;

% --- LIMITES DO ESPAÇO DE BUSCA ---
lb = 0.1;   % Área mínima de cada barra [in²]
ub = 20.0;  % Área máxima de cada barra [in²]

% --- CONFIG DO SOLVER ---
n_barras = 10;
for i = 1:n_barras
    config_vars(i).tipo = 'R';   % Área é variável contínua (Real)
    config_vars(i).min  = lb;    % Limite inferior do espaço de busca
    config_vars(i).max  = ub;    % Limite superior do espaço de busca
end

% --- PARÂMETROS DO PROBLEMA (passados para a função de custo) ---
sigma_max = 25.0;   % Tensão admissível máxima [ksi] — restrição estrutural

fprintf('============================================================\n');
fprintf(' OTIMIZAÇÃO DA TRELIÇA DE 10 BARRAS — PSO-RID\n');
fprintf(' Variáveis : %d áreas  |  Intervalo: [%.1f, %.1f] in²\n', n_barras, lb, ub);
fprintf(' Restrição : |sigma| <= %.1f ksi\n', sigma_max);
fprintf('============================================================\n\n');

% --- MULTI-START ---
n_runs             = 5;    % Número de execuções independentes
melhor_peso_global = inf;  % Melhor peso encontrado entre todos os runs
melhor_sol_global  = [];   % Melhor vetor de áreas encontrado
todos_historicos   = {};   % Histórico de convergência de cada run

% Função anônima: "congela" sigma_max para o solver só precisar passar x
funcao_custo = @(x) custo_truss(x, sigma_max);

for run = 1:n_runs
    fprintf('--- Run %d/%d ---\n', run, n_runs);

    % Chama o solver PSO-RID genérico
    % Retorna: solução (vetor de áreas), custo (peso+penalidade), histórico
    [sol, custo, hist] = pso_rid_generico(funcao_custo, config_vars);

    todos_historicos{run} = hist;   % Guarda o histórico desta run

    % Atualiza o melhor global se esta run foi melhor
    if custo < melhor_peso_global
        melhor_peso_global = custo;
        melhor_sol_global  = sol;
    end

    fprintf('  Run %d: peso = %.4f lb  |  áreas = [', run, custo);
    fprintf('%.3f ', sol);
    fprintf(']\n\n');
end

% --- VALIDAÇÃO: recalcula peso e tensões da melhor solução ---
[peso_final, sigma_final] = truss(melhor_sol_global);

fprintf('============================================================\n');
fprintf(' RESULTADO FINAL (melhor de %d runs)\n', n_runs);
fprintf('------------------------------------------------------------\n');
fprintf(' Peso total         : %.4f lb\n', peso_final);
fprintf(' Tensão admissível  : %.1f ksi\n', sigma_max);
fprintf('------------------------------------------------------------\n');
fprintf(' Barra  |  Área (in²)  |  Tensão (ksi)  |  Status\n');
fprintf('------------------------------------------------------------\n');
for i = 1:n_barras
    status = '  OK';
    if abs(sigma_final(i)) > sigma_max
        status = '  VIOLAÇÃO!';
    end
    fprintf('   %2d   |   %7.4f    |   %10.4f   | %s\n', ...
            i, melhor_sol_global(i), sigma_final(i), status);
end
fprintf('============================================================\n');

% --- PLOTS ---
figure('Name', 'PSO-RID — Treliça 10 Barras', 'NumberTitle', 'off', ...
       'Position', [100 100 1000 500]);

% Subplot 1: Convergência de todos os runs (escala log)
subplot(1,3,1);
cores = lines(n_runs);
hold on;
for run = 1:n_runs
    h = todos_historicos{run};
    semilogy(h, 'Color', [cores(run,:), 0.7], 'LineWidth', 1.2);
end
xlabel('Iteração'); ylabel('Custo (escala log)');
title('Convergência por Run');
labels_runs = arrayfun(@(x) sprintf('Run %d', x), 1:n_runs, 'UniformOutput', false);
legend(labels_runs, 'Location', 'northeast', 'FontSize', 7);
grid on; box on;

% Subplot 2: Melhor run em destaque
[~, idx_melhor_run] = min(cellfun(@(h) h(end), todos_historicos));
subplot(1,3,2);
h_melhor = todos_historicos{idx_melhor_run};
semilogy(h_melhor, 'b-', 'LineWidth', 2);
xlabel('Iteração'); ylabel('Custo (escala log)');
title(sprintf('Melhor Run (%d) — Peso: %.2f lb', idx_melhor_run, melhor_peso_global));
grid on; box on;

% Subplot 3: Áreas ótimas por barra (gráfico de barras)
subplot(1,3,3);
bar(melhor_sol_global, 'FaceColor', [0.2 0.5 0.8]);
hold on;
yline(lb, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('Mín (%.1f)', lb));
yline(ub, 'g--', 'LineWidth', 1.2, 'DisplayName', sprintf('Máx (%.1f)', ub));
xlabel('Barra'); ylabel('Área (in²)');
title('Áreas Ótimas por Barra');
legend('Área ótima', 'Limite inf.', 'Limite sup.', 'Location', 'northeast', 'FontSize', 7);
xticks(1:n_barras);
grid on; box on;

sgtitle('Otimização da Treliça de 10 Barras — PSO-RID', ...
        'FontSize', 13, 'FontWeight', 'bold');

end


% FUNÇÃO DE CUSTO DA TRELIÇA
%   J         — valor da função objetivo (peso [lb])
%   violacao  — medida de infeasibilidade (0 = solução viável)
function [J, violacao] = custo_truss(area, sigma_max)

% Chama o solver de elementos finitos da treliça
% Retorna: peso total [lb] e tensão em cada uma das 10 barras [ksi]
[weight, Sigma] = truss(area);

% --- CÁLCULO DA VIOLAÇÃO (Método de Deb) ---
% uma solução viável (violacao=0) sempre vence uma inviável (violacao>0),
violacao = 0;
for i = 1:length(Sigma)
    excesso = abs(Sigma(i)) - sigma_max;   % Quanto a tensão ultrapassa o limite
    if excesso > 0
        violacao = violacao + excesso^2;   % Penalidade quadrática: cresce rápido com a violação
    end
end

% --- FUNÇÃO OBJETIVO ---
J = weight;

end