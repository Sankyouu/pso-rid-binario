function [weight, Sigma, u_livre] = fem_truss_linear_pre_bloco2(area)
% FEM_TRUSS_LINEAR Análise Elástica Linear Clássica da Treliça de 10 Barras
% Benchmark: Hadi & Alvani (2003)
%
% Entrada:
%   area     : vetor com áreas transversais das 10 barras [mm²]
%
% Saídas:
%   weight   : massa total da estrutura [kg]
%   Sigma    : tensões axiais em cada barra [MPa]
%   u_livre  : deslocamentos nodais nos DOFs livres (DOFs 1 a 8) [mm]

% 0. Validação de Entrada
if numel(area) ~= 10
    error('fem_truss_linear:InvalidInput', ...
        'O vetor de entrada "area" deve conter exatamente 10 elementos (uma área para cada barra da treliça). Recebido: %d elementos.', numel(area));
end

% 1. Propriedades do Material (Alumínio do Benchmark Hadi 2003)
E    = 6.89e4;       % Módulo de elasticidade [MPa]
dens = 2.770e-6;     % Densidade [kg/mm³]

% 2. Geometria: Coordenadas Iniciais dos Nós [mm]
%   nó    1  2  3  4  5  6
nodes = [2  2  1  1  0  0 ;
         1  0  1  0  1  0 ] * 9144; % [mm]

% 3. Topologia: Conectividade dos 10 Elementos [nó_inicial; nó_final]
elements = [5  3  6  4  3  1  5  6  3  4 ;
            3  1  4  2  4  2  4  3  2  1];

n_el    = 10;
n_nodes = 6;
s_dof   = 12;

% 4. Condições de Contorno (Apoios fixos nos nós 5 e 6)
apoios      = [9 10 11 12];
dofs_livres = 1:8;

% 5. Forças Externas Aplicadas [N]
% 100 kip = 444.822,16 N nos nós 2 e 4 (vertical para baixo)
P_kip = 100 * 4.4482216152605 * 1000;
F = zeros(s_dof, 1);
F(4) = -P_kip; % nó 2, vertical y [N]
F(8) = -P_kip; % nó 4, vertical y [N]

% 6. Montagem da Matriz de Rigidez Global [K]
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
    c(n) = dx / L(n);
    s(n) = dy / L(n);

    ke = (area(n)*E/L(n)) * ...
        [  c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
           c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
          -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
          -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    idx = [2*i-1, 2*i, 2*j-1, 2*j];
    kk(idx, idx) = kk(idx, idx) + ke;
end

% 7. Condições de Contorno
kk_bc = kk;
F_bc  = F;
for dof = apoios
    kk_bc(dof, :)   = 0;
    kk_bc(:, dof)   = 0;
    kk_bc(dof, dof) = 1.0;
    F_bc(dof)       = 0.0;
end

% 8. Solução Linear
u = kk_bc \ F_bc;
u_livre = u(dofs_livres);

% 9. Deformações, Tensões e Massa Total
eps = zeros(n_el, 1);
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);
    idx = [2*i-1, 2*i, 2*j-1, 2*j];
    d = [-c(n), -s(n), c(n), s(n)] * u(idx);
    eps(n) = d / L(n);
end

Sigma  = E * eps;
weight = sum(area(:) .* L .* dens);

end
