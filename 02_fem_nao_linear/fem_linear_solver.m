function [weight, Sigma, u_livre, diag_out] = fem_linear_solver(problema, areas)
arguments
    problema (1,1) struct
    areas    (1,:) double {mustBePositive}
end
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
%   areas    : 1 x n_el          areas das secoes transversais. Validado pelo
%              bloco "arguments" como numerico positivo; o tamanho exato so e
%              conhecido apos ler problema.elements, entao essa checagem
%              continua no corpo da funcao.
%
% SAIDAS
%   weight   : massa total
%   Sigma    : n_el x 1  tensoes axiais (tracao > 0)
%   u_livre  : deslocamentos nos GDLs livres
%   diag_out : struct com u, L0, cossenos diretores, GDLs livres e P_axial
%
% See also fem_nao_linear_solver

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
assert(numel(problema.F_total) == s_dof, 'fem_linear_solver:forcaIncompativel', ...
    'F_total deve ter %d componentes. Recebido: %d.', s_dof, numel(problema.F_total));

F = problema.F_total(:);
dofs_livres = setdiff(1:s_dof, apoios);

% -------------------------------------------------------------------------
% 1. GEOMETRIA E MAPAS DE INDICE — pre-computados uma unica vez
% -------------------------------------------------------------------------
% Mesmo padrao vetorizado de fem_nao_linear_solver.m: monta [KE] via
% accumarray sobre indices lineares pre-computados, em vez de somar uma 4x4
% por barra dentro de um laco for. Aqui a matriz e so [KE] (sem termo
% geometrico) e a montagem acontece uma unica vez — nao ha laco de
% Newton-Raphson.

no_i = elements(1,:);
no_j = elements(2,:);

dx = nodes(1,no_j) - nodes(1,no_i);   % 1 x n_el
dy = nodes(2,no_j) - nodes(2,no_i);   % 1 x n_el
Lr = hypot(dx, dy);                    % 1 x n_el
L  = Lr.';                             % n_el x 1 (mantem o contrato de saida)

assert(all(L > 0), 'fem_linear_solver:barraDegenerada', ...
    'Ha barra(s) com comprimento nulo.');

cr = dx ./ Lr;    % 1 x n_el, cos(phi)
sr = dy ./ Lr;    % 1 x n_el, sin(phi)
c  = cr.';        % n_el x 1 (mantem o contrato de saida)
s  = sr.';        % n_el x 1

edof     = [2*no_i-1; 2*no_i; 2*no_j-1; 2*no_j];   % 4 x n_el
edof_lin = edof(:); %#ok<NASGU> — mantido por simetria com fem_nao_linear_solver.m

lin_i = repmat(edof, 4, 1);        % indice de LINHA das 16 posicoes
lin_j = kron(edof, ones(4,1));     % indice de COLUNA das 16 posicoes
linKT = lin_i(:) + (lin_j(:) - 1) * s_dof;

% -------------------------------------------------------------------------
% 2. MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL [KE]
%    (configuracao INICIAL — sem atualizacao geometrica)
% -------------------------------------------------------------------------
% A 4x4 de cada elemento tem so tres valores distintos (mesma estrutura de
% fem_nao_linear_solver.m, aqui sem o termo geometrico [KG]):
%
%     [  A   B  -A  -B ]      A = ea*cos^2
%     [  B   C  -B  -C ]      B = ea*cos*sin
%     [ -A  -B   A   B ]      C = ea*sin^2
%     [ -B  -C   B   C ]

ea = (areas .* E) ./ Lr;

A = ea .* cr.^2;
B = ea .* cr .* sr;
C = ea .* sr.^2;

Z = [ A;  B; -A; -B; ...
      B;  C; -B; -C; ...
     -A; -B;  A;  B; ...
     -B; -C;  B;  C];

kk = reshape(accumarray(linKT, Z(:), [s_dof*s_dof, 1]), s_dof, s_dof);

% -------------------------------------------------------------------------
% 3. CONDICOES DE CONTORNO
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
% 4. SOLUCAO
% -------------------------------------------------------------------------
u       = kk_bc \ F_bc;
u_livre = u(dofs_livres);

% -------------------------------------------------------------------------
% 5. DEFORMACOES, TENSOES E MASSA
% -------------------------------------------------------------------------
% Alongamento projetado sobre o eixo da barra, vetorizado para todas as
% barras de uma vez: d_n = c_n*(u_{2j-1}-u_{2i-1}) + s_n*(u_{2j}-u_{2i}).
% u e coluna (s_dof x 1); indexar por no_i/no_j (linhas) devolve coluna na
% MESMA orientacao de u — sem necessidade de transpor.

eps = (c .* (u(2*no_j-1) - u(2*no_i-1)) + ...
       s .* (u(2*no_j)   - u(2*no_i))) ./ L;

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
