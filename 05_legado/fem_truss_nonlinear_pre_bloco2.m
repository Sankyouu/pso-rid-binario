function [weight, Sigma, u_livre, diag_out] = fem_truss_nonlinear_pre_bloco2(area)
% FEM_TRUSS_NONLINEAR Análise Geometricamente Não Linear da Treliça de 10 Barras
% Benchmark: Hadi & Alvani (2003)
% Método: Newton-Raphson incremental com matrizes [KT] = [KE] + [KG]
%
% Entrada:
%   area     : vetor com áreas transversais das 10 barras [mm²]
%
% Saídas:
%   weight   : massa total da estrutura [kg]
%   Sigma    : tensões axiais em cada barra [MPa] (tração > 0, compressão < 0)
%   u_livre  : deslocamentos nodais nos DOFs livres (DOFs 1 a 8) [mm]
%   diag_out : (opcional) struct com histórico de deslocamentos, forças e nós deformados

areas_hadi = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65];

% 0. Validação de Entrada
if numel(area) ~= 10
    error('fem_truss_nonlinear:InvalidInput', ...
        'O vetor de entrada "area" deve conter exatamente 10 elementos (uma área para cada barra da treliça). Recebido: %d elementos.', numel(area));
end

% 1. Parâmetros de Controle do Newton-Raphson Incremental
n_inc    = 10;       % Número de incrementos de carga
max_iter = 50;       % Máximo de iterações por incremento
tol      = 1e-6;     % Tolerância de convergência (resíduo relativo)

% 2. Propriedades do Material (Alumínio do Benchmark Hadi 2003)
E    = 6.89e4;       % Módulo de elasticidade [MPa]
dens = 2.770e-6;     % Densidade [kg/mm³]

% 3. Geometria: Coordenadas Iniciais dos Nós [mm]
%   nó 5 ------- nó 3 ------- nó 1
%        \       |       \       |
%         \      |        \      |
%   nó 6 ------- nó 4 ------- nó 2
%      ←── 9144 mm ──→←── 9144 mm ──→
%      ↑ 9144 mm ↓
%
%   nó    1  2  3  4  5  6
nodes0 = [2  2  1  1  0  0 ;
          1  0  1  0  1  0 ] * 9144; % [mm]

% 4. Topologia: Conectividade dos 10 Elementos [nó_inicial; nó_final]
elements = [5  3  6  4  3  1  5  6  3  4 ;
            3  1  4  2  4  2  4  3  2  1];

n_el    = 10;
n_nodes = 6;
s_dof   = 12;

% 5. Condições de Contorno (Apoios fixos nos nós 5 e 6)
apoios      = [9 10 11 12];
dofs_livres = 1:8;

% 6. Forças Externas Aplicadas [N]
% Fator exato de conversão lbf -> N (NIST): 4.4482216152605 N/lbf
% 100 kip = 444.822,16 N nos nós 2 e 4 (vertical para baixo)
P_kip = 100 * 4.4482216152605 * 1000;
F_total = zeros(s_dof, 1);
F_total(4) = -P_kip; % nó 2, vertical y [N]
F_total(8) = -P_kip; % nó 4, vertical y [N]

% 7. Inicialização
u         = zeros(s_dof, 1);
P_axial   = zeros(n_el, 1);
nodes_cur = nodes0;

% Comprimentos iniciais das barras (L0)
L0 = zeros(n_el, 1);
for n = 1:n_el
    ni = elements(1,n); nj = elements(2,n);
    dx = nodes0(1,nj) - nodes0(1,ni);
    dy = nodes0(2,nj) - nodes0(2,ni);
    L0(n) = sqrt(dx^2 + dy^2);
end

% 8. Loop de Incrementos de Carga (Newton-Raphson)
dF = F_total / n_inc;

for inc = 1:n_inc
    R = dF;

    for iter = 1:max_iter
        % Montagem de KT = KE + KG
        KT = zeros(s_dof, s_dof);

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            xi = nodes_cur(1,ni); yi = nodes_cur(2,ni);
            xj = nodes_cur(1,nj); yj = nodes_cur(2,nj);

            dx = xj - xi;
            dy = yj - yi;
            Ln = sqrt(dx^2 + dy^2);

            lam = dx / Ln;
            mu  = dy / Ln;

            % KE (Rigidez elástica do elemento na config corrente)
            ke = (area(n) * E / Ln) * ...
                [  lam*lam   lam*mu  -lam*lam  -lam*mu ;
                   lam*mu    mu*mu   -lam*mu   -mu*mu  ;
                  -lam*lam  -lam*mu   lam*lam   lam*mu ;
                  -lam*mu   -mu*mu    lam*mu    mu*mu  ];

            % KG (Rigidez geométrica do elemento)
            Pn = P_axial(n);
            kg = (Pn / Ln) * ...
                [  mu*mu   -lam*mu  -mu*mu    lam*mu ;
                  -lam*mu   lam*lam  lam*mu  -lam*lam;
                  -mu*mu    lam*mu   mu*mu   -lam*mu ;
                   lam*mu  -lam*lam -lam*mu   lam*lam];

            idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
            KT(idx, idx) = KT(idx, idx) + ke + kg;
        end

        % Condições de contorno
        KT_bc = KT;
        R_bc  = R;
        for dof = apoios
            KT_bc(dof, :)   = 0;
            KT_bc(:, dof)   = 0;
            KT_bc(dof, dof) = 1.0;
            R_bc(dof)       = 0.0;
        end

        % Solução do passo N-R
        du = KT_bc \ R_bc;
        u  = u + du;

        % Atualização da geometria corrente
        for nd = 1:n_nodes
            nodes_cur(1,nd) = nodes0(1,nd) + u(2*nd-1);
            nodes_cur(2,nd) = nodes0(2,nd) + u(2*nd);
        end

        % Atualização das forças axiais (deformação de engenharia)
        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            dx0 = nodes0(1,nj) - nodes0(1,ni);
            dy0 = nodes0(2,nj) - nodes0(2,ni);
            L_ini = sqrt(dx0^2 + dy0^2);

            dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
            dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
            L_cor = sqrt(dxc^2 + dyc^2);

            eps_n = (L_cor - L_ini) / L_ini;
            P_axial(n) = area(n) * E * eps_n;
        end

        % Forças internas e resíduo
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

            idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
            felem = P_axial(n) * [-lam; -mu; lam; mu];
            F_int(idx) = F_int(idx) + felem;
        end

        R = F_ext_acc - F_int;
        R(apoios) = 0;

        norma_R = norm(R);
        norma_F = norm(F_ext_acc);
        if norma_F < 1e-12, norma_F = 1.0; end

        if (norma_R / norma_F) < tol
            break;
        end
    end
end

% 9. Pós-processamento
u_livre = u(dofs_livres);
Sigma   = P_axial ./ area(:);
weight  = sum(area(:) .* L0 .* dens);

if nargout > 3
    diag_out.u         = u;
    diag_out.P_axial   = P_axial;
    diag_out.nodes_cur = nodes_cur;
    diag_out.L0        = L0;
end

end
