function str = relatorio_comparativo(nome_caso, best_sol, best_peso, ref_sol, ref_peso, Sigma, sigma_max, u_livre, d_max, viol)
% RELATORIO_COMPARATIVO Gera relatório estruturado de resultados em texto

if nargin < 10, viol = 0; end

n_barras = length(best_sol);

fprintf('\n============================================================\n');
fprintf(' RELATÓRIO DE RESULTADOS: %s\n', upper(nome_caso));
fprintf('============================================================\n');

% Tabela de Áreas
fprintf(' ÁREAS ÓTIMAS ENCONTRADAS [mm²]:\n');
fprintf(' ------------------------------------------------------------\n');
fprintf(' Barra |   PSO-RID   |   Referência   |  Diferença\n');
fprintf(' ------------------------------------------------------------\n');
for i = 1:n_barras
    if ~isempty(ref_sol) && length(ref_sol) >= i
        fprintf('  %2d   |  %9.1f  |   %9.1f    |  %+9.1f\n', ...
                i, best_sol(i), ref_sol(i), best_sol(i) - ref_sol(i));
    else
        fprintf('  %2d   |  %9.1f  |        -       |       -\n', i, best_sol(i));
    end
end
fprintf(' ------------------------------------------------------------\n');

% Comparativo de Peso
fprintf(' Peso Total PSO-RID : %.2f kg\n', best_peso);
if ~isempty(ref_peso) && ~isnan(ref_peso)
    dif_kg  = best_peso - ref_peso;
    dif_pct = 100 * dif_kg / ref_peso;
    fprintf(' Peso Referência    : %.2f kg\n', ref_peso);
    fprintf(' Diferença          : %+0.2f kg (%+.2f%%)\n', dif_kg, dif_pct);
end
fprintf(' Violação Total     : %.6e\n', viol);
fprintf(' ------------------------------------------------------------\n');

% Tensões
if ~isempty(Sigma)
    fprintf(' TENSÕES AXIAIS NAS BARRAS [MPa] (sigma_max = %.2f MPa):\n', sigma_max);
    fprintf(' Barra |   Tensão [MPa]  | |sigma|/sigma_max | Status\n');
    fprintf(' ------------------------------------------------------------\n');
    for i = 1:n_barras
        taxa = abs(Sigma(i)) / sigma_max;
        status = 'OK';
        if taxa > 1.0, status = 'VIOLOU'; end
        fprintf('  %2d   |   %+11.3f   |      %6.3f       | %s\n', ...
                i, Sigma(i), taxa, status);
    end
    fprintf(' ------------------------------------------------------------\n');
end

% Deslocamentos
if ~isempty(u_livre) && nargin >= 9 && ~isempty(d_max)
    fprintf(' DESLOCAMENTOS DOS DOFs LIVRES [mm] (d_max = %.2f mm):\n', d_max);
    fprintf(' DOF   | Deslocamento [mm] |   |u|/d_max   | Status\n');
    fprintf(' ------------------------------------------------------------\n');
    for j = 1:length(u_livre)
        taxa = abs(u_livre(j)) / d_max;
        status = 'OK';
        if taxa > 1.0, status = 'VIOLOU'; end
        fprintf('  %2d   |    %+11.4f    |    %6.3f     | %s\n', ...
                j, u_livre(j), taxa, status);
    end
    fprintf(' Deslocamento máximo absoluto: %.4f mm\n', max(abs(u_livre)));
end

fprintf('============================================================\n\n');

end
