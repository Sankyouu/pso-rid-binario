function [weight, Sigma, u_livre, diag_out] = fem_nao_linear_solver(problema, areas, opcoes)
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

% Parametros de controle do Newton-Raphson
if nargin < 3 || isempty(opcoes), opcoes = struct(); end
n_inc    = get_opcao(opcoes, 'n_inc',      10);
max_iter = get_opcao(opcoes, 'max_iter',   50);
tol      = get_opcao(opcoes, 'tol',      1e-6);

dofs_livres = setdiff(1:s_dof, apoios);

% -------------------------------------------------------------------------
% 1. COMPRIMENTOS INICIAIS (L0)
% -------------------------------------------------------------------------

L0 = zeros(n_el, 1);
for n = 1:n_el
    ni = elements(1,n);  nj = elements(2,n);
    dx = nodes0(1,nj) - nodes0(1,ni);
    dy = nodes0(2,nj) - nodes0(2,ni);
    L0(n) = sqrt(dx^2 + dy^2);
end
assert(all(L0 > 0), 'fem_nao_linear_solver:barraDegenerada', ...
    'Ha barra(s) com comprimento inicial nulo.');

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

        KT = zeros(s_dof, s_dof);

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            dx = nodes_cur(1,nj) - nodes_cur(1,ni);
            dy = nodes_cur(2,nj) - nodes_cur(2,ni);
            Ln = sqrt(dx^2 + dy^2);

            % Cossenos diretores na configuracao corrente
            % [HA2003] Eq. (14): lambda = cos(phi), mu = sin(phi)

            lam = dx / Ln;
            mu  = dy / Ln;

            % --- [KE] rigidez elastica (mesma forma da analise linear,
            %     porem com a geometria corrente) ---

            ke = (areas(n) * E / Ln) * ...
                [  lam*lam   lam*mu  -lam*lam  -lam*mu ;
                   lam*mu    mu*mu   -lam*mu   -mu*mu  ;
                  -lam*lam  -lam*mu   lam*lam   lam*mu ;
                  -lam*mu   -mu*mu    lam*mu    mu*mu  ];

            % --- [KG] rigidez geometrica ---
            % Eq. (14) [HA2003]:
            %   KG = (P/L) * [  mu^2   -lam*mu  -mu^2    lam*mu ;
            %                  -lam*mu  lam^2    lam*mu -lam^2  ;
            %                  -mu^2    lam*mu   mu^2   -lam*mu ;
            %                   lam*mu -lam^2   -lam*mu  lam^2  ]
            % onde P e a carga axial intermediaria do elemento no estagio
            % de carregamento corrente e L o comprimento do elemento.

            Pn = P_axial(n);
            kg = (Pn / Ln) * ...
                [  mu*mu   -lam*mu  -mu*mu    lam*mu ;
                  -lam*mu   lam*lam  lam*mu  -lam*lam;
                  -mu*mu    lam*mu   mu*mu   -lam*mu ;
                   lam*mu  -lam*lam -lam*mu   lam*lam];

            idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
            KT(idx, idx) = KT(idx, idx) + ke + kg;   % Eq. (13)
        end

        % --- Condicoes de contorno: zera linha/coluna e poe 1 na diagonal

        KT_bc = KT;
        R_bc  = R;
        for dof = apoios
            KT_bc(dof, :)   = 0;
            KT_bc(:, dof)   = 0;
            KT_bc(dof, dof) = 1.0;
            R_bc(dof)       = 0.0;
        end

        % --- Passo de Newton-Raphson: Eq. (18) [HA2003]
        %     [Kt]_i {du} = {R}_i

        du = KT_bc \ R_bc;
        u  = u + du;

        % --- Atualizacao da configuracao corrente

        for nd = 1:n_nodes
            nodes_cur(1,nd) = nodes0(1,nd) + u(2*nd-1);
            nodes_cur(2,nd) = nodes0(2,nd) + u(2*nd);
        end

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

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);
            dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
            dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
            L_cor = sqrt(dxc^2 + dyc^2);

            eps_n      = (L_cor - L0(n)) / L0(n);
            P_axial(n) = areas(n) * E * eps_n;
        end

        % --- Forcas internas na configuracao corrente

        F_ext_acc = inc * dF;
        F_int     = zeros(s_dof, 1);

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);
            dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
            dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
            Ln  = sqrt(dxc^2 + dyc^2);
            lam = dxc / Ln;
            mu  = dyc / Ln;

            % Equilibrio do elemento: q = P * [-lam; -mu; +lam; +mu]

            idx   = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
            felem = P_axial(n) * [-lam; -mu; lam; mu];
            F_int(idx) = F_int(idx) + felem;
        end

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
