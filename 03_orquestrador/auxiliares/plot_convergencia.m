function fig = plot_convergencia(todos_historicos, ref_peso, titulo_str)
% PLOT_CONVERGENCIA Gera visualizações padronizadas de convergência para múltiplos runs
%
% Entradas:
%   todos_historicos : cell array contendo o histórico de cada run (vetores coluna)
%   ref_peso         : valor de referência do peso (escalar) para traçar linha guia
%   titulo_str       : título principal da figura

if nargin < 3 || isempty(titulo_str)
    titulo_str = 'Curvas de Convergência PSO-RID';
end

n_runs = length(todos_historicos);
fig = figure('Name', titulo_str, 'NumberTitle', 'off', 'Position', [150 150 950 480]);

% Subplot 1: Todos os runs
subplot(1, 2, 1);
cores = lines(n_runs);
hold on;
for r = 1:n_runs
    semilogy(todos_historicos{r}, 'Color', [cores(r,:), 0.6], 'LineWidth', 1.2);
end

if nargin >= 2 && ~isempty(ref_peso) && ~isnan(ref_peso)
    yline(ref_peso, 'k--', 'LineWidth', 1.5);
    legend_entries = arrayfun(@(x) sprintf('Run %d', x), 1:n_runs, 'UniformOutput', false);
    legend_entries{end+1} = 'Ref. Literatura';
    legend(legend_entries, 'Location', 'northeast', 'FontSize', 7);
end

xlabel('Iteração', 'FontWeight', 'bold');
ylabel('Peso [kg] (escala log)', 'FontWeight', 'bold');
title('Convergência de Todos os Runs', 'FontWeight', 'bold');
grid on; box on;

% Subplot 2: Melhor Run
custos_finais = cellfun(@(h) h(end), todos_historicos);
[melhor_peso, idx_melhor] = min(custos_finais);

subplot(1, 2, 2);
semilogy(todos_historicos{idx_melhor}, 'b-', 'LineWidth', 2.0);
hold on;

if nargin >= 2 && ~isempty(ref_peso) && ~isnan(ref_peso)
    yline(ref_peso, 'r--', 'LineWidth', 1.5);
    legend(sprintf('Melhor (Run %d: %.2f kg)', idx_melhor, melhor_peso), 'Ref. Literatura', 'Location', 'northeast');
else
    legend(sprintf('Melhor (Run %d: %.2f kg)', idx_melhor, melhor_peso), 'Location', 'northeast');
end

xlabel('Iteração', 'FontWeight', 'bold');
ylabel('Peso [kg] (escala log)', 'FontWeight', 'bold');
title(sprintf('Melhor Execução (Run %d)', idx_melhor), 'FontWeight', 'bold');
grid on; box on;

sgtitle(titulo_str, 'FontSize', 12, 'FontWeight', 'bold');

end
