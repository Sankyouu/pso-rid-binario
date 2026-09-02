function [weight, Sigma, u_livre, diag_out] = fem_linear_solver(problema, areas)
% FEM_LINEAR_SOLVER  Analise elastica linear de trelicas planas.
%
% BLOCO 2a do projeto CPIO III — SOLVER GENERICO (contraparte linear).
%
% Como o solver nao linear, este arquivo NAO conhece nenhum problema
% especifico: tudo chega pelo struct 'problema'.
%
% -------------------------------------------------------------------------
% REFERENCIA
% -------------------------------------------------------------------------
% [HA2003] Hadi & Alvani (2003), Paper 37, Civil-Comp Press.
%          Sec. 4: "During loading, initially no change in geometry is
%          registered, KG = 0, rendering this stage similar to a linear
%          analysis."
%
% Ou seja, a analise linear equivale a montar apenas [KE] na configuracao
% INICIAL e resolver [KE]{u} = {F} em um unico passo. Esta equivalencia e
% verificada por teste unitario: para carga tendendo a zero, o solver nao
% linear deve convergir para este resultado.
%
% -------------------------------------------------------------------------
% ENTRADAS
% -------------------------------------------------------------------------
%   problema : struct (mesmo formato do solver nao linear)
%       .nodes0    2 x n_nodes   coordenadas [x; y]
%       .elements  2 x n_el      conectividade [no_inicial; no_final]
%       .apoios    1 x n_ap      GDLs restringidos
%       .F_total   s_dof x 1     forcas externas
%       .E         escalar       modulo de elasticidade
%       .dens      escalar       densidade
%   areas    : 1 x n_el          areas das secoes transversais
%
% SAIDAS
%   weight   : massa total
%   Sigma    : n_el x 1  tensoes axiais (tracao > 0)
%   u_livre  : deslocamentos nos GDLs livres
%   diag_out : struct com u, L0, cossenos diretores, GDLs livres e P_axial

% -------------------------------------------------------------------------
% 0. VALIDACAO
% -------------------------------------------------------------------------
campos_obrig = {'nodes0', 'elements', 'apoios', 'F_total', 'E', 'dens'};
for i = 1:numel(campos_obrig)
    assert(isfield(problema, campos_obrig{i}), ...
        'fem_linear_solver:campoAusente', ...
        'O struct "problema" deve conter o campo "%s".', campos_obrig{i});
end

nodes    = problema.nodes0;
elements = problema.elements;
apoios   = problema.apoios(:)';
E        = problema.E;
dens     = problema.dens;

n_nodes = size(nodes, 2);
n_el    = size(elements, 2);
s_dof   = 2 * n_nodes;

assert(size(nodes, 1) == 2, 'fem_linear_solver:geometria2D', ...
    'nodes0 deve ter 2 linhas (trelica plana): [x; y].');
assert(numel(areas) == n_el, 'fem_linear_solver:areasIncompativeis', ...
    'O vetor "areas" deve ter %d elementos. Recebido: %d.', n_el, numel(areas));
assert(all(areas(:) > 0), 'fem_linear_solver:areaNaoPositiva', ...
    'Todas as areas devem ser estritamente positivas.');
assert(numel(problema.F_total) == s_dof, 'fem_linear_solver:forcaIncompativel', ...
    'F_total deve ter %d componentes. Recebido: %d.', s_dof, numel(problema.F_total));

F = problema.F_total(:);
dofs_livres = setdiff(1:s_dof, apoios);

% -------------------------------------------------------------------------
% 1. MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL [KE]
%    (configuracao INICIAL — sem atualizacao geometrica)
% -------------------------------------------------------------------------
kk = zeros(s_dof, s_dof);
L  = zeros(n_el, 1);
c  = zeros(n_el, 1);
s  = zeros(n_el, 1);

for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);

    dx = nodes(1,j) - nodes(1,i);
    dy = nodes(2,j) - nodes(2,i);
    L(n) = sqrt(dx^2 + dy^2);
    assert(L(n) > 0, 'fem_linear_solver:barraDegenerada', ...
        'Barra %d tem comprimento nulo.', n);

    c(n) = dx / L(n);      % cos(phi)
    s(n) = dy / L(n);      % sin(phi)

    ke = (areas(n)*E/L(n)) * ...
        [  c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
           c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
          -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
          -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    idx = [2*i-1, 2*i, 2*j-1, 2*j];
    kk(idx, idx) = kk(idx, idx) + ke;
end

% -------------------------------------------------------------------------
% 2. CONDICOES DE CONTORNO
% -------------------------------------------------------------------------
kk_bc = kk;
F_bc  = F;
for dof = apoios
    kk_bc(dof, :)   = 0;
    kk_bc(:, dof)   = 0;
    kk_bc(dof, dof) = 1.0;
    F_bc(dof)       = 0.0;
end

% -------------------------------------------------------------------------
% 3. SOLUCAO
% -------------------------------------------------------------------------
u       = kk_bc \ F_bc;
u_livre = u(dofs_livres);

% -------------------------------------------------------------------------
% 4. DEFORMACOES, TENSOES E MASSA
% -------------------------------------------------------------------------
eps = zeros(n_el, 1);
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);
    idx = [2*i-1, 2*i, 2*j-1, 2*j];
    % Alongamento projetado sobre o eixo da barra
    d = [-c(n), -s(n), c(n), s(n)] * u(idx);
    eps(n) = d / L(n);
end

Sigma  = E * eps;
weight = sum(areas(:) .* L .* dens);

if nargout > 3
    diag_out.u           = u;
    diag_out.L0          = L;
    diag_out.cos_dir     = c;
    diag_out.sin_dir     = s;
    diag_out.dofs_livres = dofs_livres;
    diag_out.P_axial     = Sigma .* areas(:);
    diag_out.K           = kk;
    diag_out.n_el        = n_el;
    diag_out.n_nodes     = n_nodes;
end

end
