% Script de Diagnóstico e Verificação dos Solvers FEM (Linear e Não Linear)
clear; clc;

% Adicionar caminhos do projeto
root_dir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(root_dir);
setup_paths(false);

% 1. VERIFICAÇÃO COM A SOLUÇÃO DE REFERÊNCIA DO BENCHMARK HADI (2003)
disp('======================================================');
disp('1. VERIFICAÇÃO COM A SOLUÇÃO DE REFERÊNCIA (HADI 2003)');
disp('======================================================');

caso = catalogo_hadi();
areas_ref = caso.ref_areas;
disp(['Peso de Referência (Hadi): ', num2str(caso.ref_peso), ' kg']);

% Chamada do solver não linear autocontido
[weight_nl, Sigma_nl, u_livre_nl, diag_nl] = fem_truss_nonlinear(areas_ref);

disp(['Peso Calculado (NL): ', num2str(weight_nl), ' kg']);
fprintf('Verificação de Peso: %.2f kg (Esperado: %.2f kg)\n', weight_nl, caso.ref_peso);

disp(' ');
disp('Deslocamentos nos DOFs Livres (NL):');
for i = 1:length(u_livre_nl)
    fprintf('DOF %d: %9.4f mm\n', i, u_livre_nl(i));
end
max_disp = max(abs(u_livre_nl));
fprintf('Deslocamento Máximo: %.4f mm (Limite: %.2f mm)\n', max_disp, caso.d_max);

disp(' ');
disp('Tensões nas Barras (NL - MPa):');
for i = 1:length(Sigma_nl)
    fprintf('Barra %2d: %9.4f MPa\n', i, Sigma_nl(i));
end
max_sigma = max(abs(Sigma_nl));
fprintf('Tensão Máxima: %.4f MPa (Limite: %.2f MPa)\n', max_sigma, caso.sigma_max);

% 2. COMPARAÇÃO LINEAR vs NÃO LINEAR
disp(' ');
disp('======================================================');
disp('2. COMPARAÇÃO LINEAR vs NÃO LINEAR');
disp('======================================================');
[weight_l, Sigma_l, u_livre_l] = fem_truss_linear(areas_ref);

disp('Deslocamentos Lineares vs Não Lineares:');
for i = 1:length(u_livre_l)
    diff_val = abs(u_livre_nl(i) - u_livre_l(i));
    pct_val = diff_val / abs(u_livre_l(i)) * 100;
    if abs(u_livre_l(i)) < 1e-10, pct_val = 0; end
    fprintf('DOF %2d: Linear=%9.4f mm, NL=%9.4f mm, Dif=%7.4f %%\n', i, u_livre_l(i), u_livre_nl(i), pct_val);
end

% 3. VERIFICAÇÃO DE EQUILÍBRIO (F_ext - F_int)
disp(' ');
disp('======================================================');
disp('3. VERIFICAÇÃO DE EQUILÍBRIO GLOBAL');
disp('======================================================');
P_kip = 100 * 4.4482216152605 * 1000;
F_ext = zeros(12, 1);
F_ext(4) = -P_kip;
F_ext(8) = -P_kip;

elements = [5  3  6  4  3  1  5  6  3  4 ;
            3  1  4  2  4  2  4  3  2  1];

u_full = diag_nl.u;
nodes_cur = diag_nl.nodes_cur;
F_int = zeros(12, 1);

for k = 1:size(elements, 2)
    ni = elements(1,k);
    nj = elements(2,k);
    dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
    dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
    Ln = sqrt(dxc^2 + dyc^2);
    lam = dxc / Ln;
    mu = dyc / Ln;
    idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
    felem = diag_nl.P_axial(k) * [-lam; -mu; lam; mu];
    F_int(idx) = F_int(idx) + felem;
end

apoios = [9 10 11 12];
Residual = F_ext - F_int;
Residual(apoios) = 0;
disp('Resíduo de Forças nos DOFs Livres:');
disp(Residual(1:8));
fprintf('Norma do Resíduo: %g N\n', norm(Residual));
