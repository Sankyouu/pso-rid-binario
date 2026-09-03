function [weight, Sigma, u_livre, diag_out] = fem_nao_linear_solver(problema, areas, opcoes)
arguments
    problema (1,1) struct
    areas    (1,:) double
    opcoes   (1,1) struct = struct()
end
% FEM_NAO_LINEAR_SOLVER  Analise geometricamente nao linear de trelicas planas.
%
% Este arquivo nao conhece nenhum problema especifico: toda a geometria,
% material, carregamento e condicoes de contorno chegam pelo struct
% 'problema' (ver 02_fem_nao_linear/problemas/ para os casos concretos).
%
% -------------------------------------------------------------------------
% REFERENCIAS TEORICAS
% -------------------------------------------------------------------------
%
% [HA2003]  Hadi, M.N.S.; Alvani, K.S. "Discrete Optimum Design of
%           Geometrically Non-Linear Trusses using Genetic Algorithms".
%           Proc. 7th Int. Conf. on the Application of Artificial
%           Intelligence to Civil and Structural Engineering, Civil-Comp
%           Press, 2003, Paper 37.
%           -> 00_docs/artigos/Discrete Optimum Design of Geometrically
%              Non-Linear Trusses using Genetic Algorithms.PDF
%
% [CRI1991] Crisfield, M.A. "Non-Linear Finite Element Analysis of Solids
%           and Structures", Vol. 1 — Essentials. John Wiley & Sons, 1991.
%           (referencia [16] de [HA2003])
%           -> 00_docs/livros/
%
% Mapeamento equacao -> local no codigo:
%   Eq. (13) [HA2003]  [KT] = [KE] + [KG] ............ secao 4 (montagem)
%   Eq. (14) [HA2003]  matriz de rigidez geometrica .. secao 4 (bloco kg)
%   Eq. (17) [HA2003]  [Kt]{d} = {F} ................. secao 4 (solucao)
%   Eq. (18) [HA2003]  [Kt]_i {du} = {dF}_i .......... secao 4 (passo N-R)
%
% [HA2003] Sec. 4.3: "In the incremental load methods, the Newton-Raphson
% method is employed on the basis of iterative procedures... An increment
% load dF_i of the total external load is applied to the structure, then the
% incremental displacements {du} is calculated by making use of the system
% tangent stiffness matrix [Kt]_i for the load increment i".
%
% -------------------------------------------------------------------------
% ENTRADAS
% -------------------------------------------------------------------------
%
%   problema : struct com os parametros do caso (ver problemas/*.m)
%       .nodes0    2 x n_nodes   coordenadas iniciais [x; y]
%       .elements  2 x n_el      conectividade [no_inicial; no_final]
%       .apoios    1 x n_ap      indices dos GDLs restringidos
%       .F_total   s_dof x 1     vetor de forcas externas totais
%       .E         escalar       modulo de elasticidade
%       .dens      escalar       densidade (massa por volume)
%
%   areas    : 1 x n_el  areas das secoes transversais
%
%   opcoes   : (opcional) struct de controle do Newton-Raphson
%       .n_inc     numero de incrementos de carga        (padrao 10)
%       .max_iter  maximo de iteracoes N-R por incremento (padrao 50)
%       .tol       tolerancia do residuo relativo         (padrao 1e-6)
%
% CONSISTENCIA DE UNIDADES: o solver e agnostico a unidades, mas elas devem
% ser coerentes entre si. Nos casos deste projeto usa-se
%   comprimento [mm] | forca [N] | tensao [MPa = N/mm^2] | massa [kg]
% o que exige densidade em [kg/mm^3].
%
% -------------------------------------------------------------------------
% SAIDAS
% -------------------------------------------------------------------------
%
%   weight   : massa total da estrutura (usa os comprimentos INICIAIS L0,
%              pois o peso e propriedade do projeto, nao do estado deformado)
%   Sigma    : n_el x 1  tensoes axiais (tracao > 0, compressao < 0)
%   u_livre  : deslocamentos nos GDLs livres
%   diag_out : struct de diagnostico
%       .u              deslocamentos em TODOS os GDLs
%       .P_axial        forcas axiais nas barras
%       .nodes_cur      coordenadas deformadas
%       .L0             comprimentos iniciais
%       .dofs_livres    indices dos GDLs livres
%       .convergiu      true se TODOS os incrementos convergiram
%       .iters_por_inc  iteracoes N-R gastas em cada incremento
%       .residuo_final  norma relativa do residuo ao fim do ultimo incremento
%
% See also fem_linear_solver

% -------------------------------------------------------------------------
% 0. VALIDACAO DE ENTRADA
% -------------------------------------------------------------------------

campos_obrig = {'nodes0', 'elements', 'apoios', 'F_total', 'E', 'dens'};
for i = 1:numel(campos_obrig)
    assert(isfield(problema, campos_obrig{i}), ...
        'fem_nao_linear_solver:campoAusente', ...
        'O struct "problema" deve conter o campo "%s".', campos_obrig{i});
end

nodes0   = problema.nodes0;
elements = problema.elements;
apoios   = problema.apoios(:)';
E        = problema.E;
dens     = problema.dens;

n_nodes = size(nodes0, 2);
n_el    = size(elements, 2);
s_dof   = 2 * n_nodes;                 % trelica plana: 2 GDLs por no

assert(size(nodes0, 1) == 2, 'fem_nao_linear_solver:geometria2D', ...
    'nodes0 deve ter 2 linhas (trelica plana): [x; y].');
assert(size(elements, 1) == 2, 'fem_nao_linear_solver:conectividade', ...
    'elements deve ter 2 linhas: [no_inicial; no_final].');
assert(numel(areas) == n_el, 'fem_nao_linear_solver:areasIncompativeis', ...
    'O vetor "areas" deve ter %d elementos (um por barra). Recebido: %d.', ...
    n_el, numel(areas));
assert(all(areas(:) > 0), 'fem_nao_linear_solver:areaNaoPositiva', ...
    'Todas as areas devem ser estritamente positivas.');
assert(numel(problema.F_total) == s_dof, 'fem_nao_linear_solver:forcaIncompativel', ...
    'F_total deve ter %d componentes (2 por no). Recebido: %d.', ...
    s_dof, numel(problema.F_total));
assert(all(apoios >= 1 & apoios <= s_dof), 'fem_nao_linear_solver:apoioInvalido', ...
    'Indices de apoio fora da faixa [1, %d].', s_dof);

F_total = problema.F_total(:);

% Parametros de controle do Newton-Raphson. O default de "opcoes" (struct
% vazio) ja vem do bloco "arguments"; get_opcao le cada campo com seu proprio
% default, porque um struct generico nao tem seus subcampos declaraveis em
% "arguments" sem trocar a convencao de chamada (opcoes so aceita struct
% pronto posicional, nao pares Nome-Valor — ver nota em fem_linear_solver.m).
n_inc    = get_opcao(opcoes, 'n_inc',      10);
max_iter = get_opcao(opcoes, 'max_iter',   50);
tol      = get_opcao(opcoes, 'tol',      1e-6);

% Mascara logica em vez de setdiff: mesmo conjunto, sem ordenar nem testar
% unicidade. setdiff custava 2% do tempo do solver (profile de 2026-09-02),
% por ser chamado uma vez a cada avaliacao da funcao objetivo.
mask_livres         = true(1, s_dof);
mask_livres(apoios) = false;
dofs_livres         = find(mask_livres);

% -------------------------------------------------------------------------
% 1. COMPRIMENTOS INICIAIS (L0)
% -------------------------------------------------------------------------

no_i = elements(1,:);
no_j = elements(2,:);
L0   = hypot(nodes0(1,no_j) - nodes0(1,no_i), ...
             nodes0(2,no_j) - nodes0(2,no_i)).';

L0r       = L0.';            % versoes LINHA, usadas nas contas vetorizadas
areas_lin = areas(:).';
assert(all(L0 > 0), 'fem_nao_linear_solver:barraDegenerada', ...
    'Ha barra(s) com comprimento inicial nulo.');

% -------------------------------------------------------------------------
% 1b. MAPAS DE INDICE — pre-computados uma unica vez
% -------------------------------------------------------------------------
% A topologia da trelica nao muda durante a analise, entao os indices de
% montagem sao constantes. Pre-computa-los tira toda a indexacao de dentro do
% laco de Newton-Raphson, que e onde estava quase metade do custo do solver.
%
%   edof_lin : GDLs de cada elemento, empilhados (4*n_el x 1)
%   linKT    : indices LINEARES das 16 posicoes da 4x4 de cada elemento
%              dentro de KT, em ordem column-major (a linha varia mais rapido)

edof     = [2*no_i-1; 2*no_i; 2*no_j-1; 2*no_j];   % 4 x n_el
edof_lin = edof(:);

lin_i = repmat(edof, 4, 1);        % indice de LINHA  das 16 posicoes
lin_j = kron(edof, ones(4,1));     % indice de COLUNA das 16 posicoes
linKT = lin_i(:) + (lin_j(:) - 1) * s_dof;

% -------------------------------------------------------------------------
% 2. INICIALIZACAO
% -------------------------------------------------------------------------

u          = zeros(s_dof, 1);   % deslocamentos acumulados
P_axial    = zeros(n_el, 1);    % forcas axiais acumuladas
nodes_cur  = nodes0;            % configuracao corrente

iters_por_inc = zeros(n_inc, 1);
convergiu_inc = false(n_inc, 1);
residuo_final = NaN;

% -------------------------------------------------------------------------
% 3. LOOP DE INCREMENTOS DE CARGA
%    [HA2003] Sec. 4.3 / Eq. (18)
% -------------------------------------------------------------------------

dF = F_total / n_inc;

for inc = 1:n_inc

    % No inicio do incremento o residuo e o proprio incremento de carga
    % (o passo anterior ja estava equilibrado).
    R = dF;

    for iter = 1:max_iter

        % -----------------------------------------------------------------
        % 4. MONTAGEM DE [KT] = [KE] + [KG]     Eq. (13) [HA2003]
        %    Avaliada na configuracao CORRENTE (nodes_cur) e com as forcas
        %    axiais correntes (P_axial).
        % -----------------------------------------------------------------

        % Geometria corrente de TODOS os elementos, de uma vez.
        dx  = nodes_cur(1,no_j) - nodes_cur(1,no_i);
        dy  = nodes_cur(2,no_j) - nodes_cur(2,no_i);
        Ln  = hypot(dx, dy);
        lam = dx ./ Ln;              % [HA2003] Eq. (14): lambda = cos(phi)
        mu  = dy ./ Ln;              %                    mu     = sin(phi)

        ll = lam .* lam;
        lm = lam .* mu;
        mm = mu  .* mu;

        ea = (areas_lin * E) ./ Ln;  % coeficiente de [KE]
        pl = P_axial.'     ./ Ln;    % coeficiente de [KG]

        % Somando [KE] e [KG] ANALITICAMENTE, a 4x4 do elemento (Eq. 13 com
        % Eq. 14) tem apenas tres valores distintos:
        %
        %     [  A   B  -A  -B ]      A = ea*lam^2 + pl*mu^2
        %     [  B   C  -B  -C ]      B = (ea - pl)*lam*mu
        %     [ -A  -B   A   B ]      C = ea*mu^2  + pl*lam^2
        %     [ -B  -C   B   C ]
        %
        % Assim nao e preciso montar e somar duas matrizes 4x4 por elemento a
        % cada iteracao de Newton-Raphson.

        A = ea.*ll + pl.*mm;
        B = (ea - pl) .* lm;
        C = ea.*mm + pl.*ll;

        % As 16 posicoes da 4x4 em ordem column-major, uma coluna por
        % elemento, acumuladas em KT pelos indices lineares pre-computados.
        Z = [ A;  B; -A; -B; ...
              B;  C; -B; -C; ...
             -A; -B;  A;  B; ...
             -B; -C;  B;  C];

        KT = reshape(accumarray(linKT, Z(:), [s_dof*s_dof, 1]), s_dof, s_dof);

        % --- Condicoes de contorno e passo de Newton-Raphson
        %     Eq. (18) [HA2003]:  [Kt]_i {du} = {R}_i
        %
        %     Resolve-se so no bloco LIVRE, o que e algebricamente identico a
        %     zerar linha e coluna de cada apoio com 1 na diagonal — e mais
        %     barato: nao copia KT nem fatora as equacoes restringidas.

        du = zeros(s_dof, 1);
        du(dofs_livres) = KT(dofs_livres, dofs_livres) \ R(dofs_livres);
        u  = u + du;

        % --- Atualizacao da configuracao corrente

        nodes_cur = nodes0 + reshape(u, 2, n_nodes);

        % --- Atualizacao das forcas axiais
        %     Medida de deformacao: DEFORMACAO DE ENGENHARIA
        %         eps = (L_corrente - L_inicial) / L_inicial
        %     N = E * A * eps
        %
        %     NOTA IMPORTANTE PARA VALIDACAO: esta escolha de medida de
        %     deformacao precisa ser levada em conta ao comparar com
        %     solucoes analiticas da literatura — formulas classicas de
        %     trelica abatida (p.ex. em [CRI1991]) costumam ser derivadas
        %     com deformacao de Green, que NAO coincide com a de engenharia
        %     fora do regime de deformacoes infinitesimais.

        dxc   = nodes_cur(1,no_j) - nodes_cur(1,no_i);
        dyc   = nodes_cur(2,no_j) - nodes_cur(2,no_i);
        L_cor = hypot(dxc, dyc);

        P_axial = (areas_lin * E .* (L_cor - L0r) ./ L0r).';

        % --- Forcas internas na configuracao corrente
        %     Equilibrio do elemento: q = P * [-lam; -mu; +lam; +mu]
        %     Reaproveita dxc/dyc/L_cor calculados logo acima.

        F_ext_acc = inc * dF;

        lamc = dxc ./ L_cor;
        muc  = dyc ./ L_cor;
        Pl   = P_axial.';

        Vf = [-Pl.*lamc; -Pl.*muc; Pl.*lamc; Pl.*muc];

        F_int = accumarray(edof_lin, Vf(:), [s_dof, 1]);

        % --- Residuo e criterio de convergencia

        R = F_ext_acc - F_int;
        R(apoios) = 0;

        norma_R = norm(R);
        norma_F = norm(F_ext_acc);
        if norma_F < 1e-12, norma_F = 1.0; end
        residuo_rel = norma_R / norma_F;

        iters_por_inc(inc) = iter;

        if residuo_rel < tol
            convergiu_inc(inc) = true;
            break;
        end
    end

    residuo_final = residuo_rel;
end

% -------------------------------------------------------------------------
% 5. POS-PROCESSAMENTO
% -------------------------------------------------------------------------

u_livre = u(dofs_livres);
Sigma   = P_axial ./ areas(:);
weight  = sum(areas(:) .* L0 .* dens);

if nargout > 3
    diag_out.u             = u;
    diag_out.P_axial       = P_axial;
    diag_out.nodes_cur     = nodes_cur;
    diag_out.L0            = L0;
    diag_out.dofs_livres   = dofs_livres;
    diag_out.convergiu     = all(convergiu_inc);
    diag_out.iters_por_inc = iters_por_inc;
    diag_out.residuo_final = residuo_final;
    diag_out.n_el          = n_el;
    diag_out.n_nodes       = n_nodes;
end

end


function val = get_opcao(s, campo, padrao)
if isfield(s, campo) && ~isempty(s.(campo))
    val = s.(campo);
else
    val = padrao;
end
end
