function [weight, Sigma, u_livre] = truss(area)

% [weight, Sigma, u_livre] = truss([19355, 65, 16129, 7742, 65, 65, 5161,
% 16129, 12903, 65]) - Caso 2 Discreto
% [weight, Sigma, u_livre] = truss([19652, 65, 13329, 11690, 65, 65, 5050, 13643, 13458, 65]) - Caso Contínuo



%   Comprimento : mm
%   Força       : N
%   Tensão      : MPa (= N/mm²)
%   Massa       : kg
%
% Referência de comparação: Tabela 1, Case 2 (catálogo 2)
%   Peso ótimo discreto = 2325.2 kg
%
% SAÍDAS:
%   weight   — Peso total [kg]
%   Sigma    — Tensão em cada barra [MPa]
%   u_livre  — Deslocamentos nodais livres [mm]  (DOFs não restringidos)

% --- PROPRIEDADES DO MATERIAL ---
E    = 6.89e4;      % Módulo de elasticidade [MPa = N/m²]
dens = 2.770e-6;    % Densidade [kg/m³]

% --- GEOMETRIA: POSIÇÃO DOS NÓS [mm] ---
%
%   5 ------- 3 ------- 1
%   |    \    |    \    |
%   6 ------- 4 ------- 2
%
%       9144mm    9144mm     (vão horizontal)
%       9144mm               (altura)
%
%   nó   1     2     3     4     5     6
nodes = [2     2     1     1     0     0  ;   % coord x [×9144 mm]
         1     0     1     0     1     0  ] * 9144; % coord y

% --- TOPOLOGIA: CONECTIVIDADE ---
%       el  1  2  3  4  5  6  7  8  9 10
elements = [5  3  6  4  3  1  5  6  3  4 ;   % nó inicial
            3  1  4  2  4  2  4  3  2  1];   % nó final

% --- DIMENSÕES DO SISTEMA ---
n_el    = size(elements, 2); % 10 elementos
ndof    = 2;                 % 2 DOFs por nó (x, y)
n_nodes = size(nodes, 2);    % 6 nós
s_dof   = n_nodes * ndof;    % 12 DOFs totais

% --- MATRIZ DE RIGIDEZ GLOBAL ---
kk = zeros(s_dof, s_dof);

% --- APOIOS ---
% Nós 5 e 6 fixos (DOFs 9,10,11,12 restritos)
apoio = [9 10 11 12];

% DOFs livres (para extrair deslocamentos relevantes)
todos_dofs = 1:s_dof;
dofs_livres = setdiff(todos_dofs, apoio);

% --- FORÇAS EXTERNAS [N] ---
% 445.4 kN nos nós 2 e 4 (verticais, para baixo)
% DOF 4  = vertical do nó 2 → -445400 N
% DOF 8  = vertical do nó 4 → -445400 N
F = [0; 0; 0; -445400; 0; 0; 0; -445400; 0; 0; 0; 0];

% MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);

    L(n) = sqrt((nodes(1,j)-nodes(1,i))^2 + (nodes(2,j)-nodes(2,i))^2);
    c(n) = (nodes(1,j) - nodes(1,i)) / L(n);
    s(n) = (nodes(2,j) - nodes(2,i)) / L(n);

    k(:,:,n) = (area(n)*E/L(n)) * ...
        [ c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
          c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
         -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
         -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    index(1) = 2*elements(1,n) - 1;
    index(2) = 2*elements(1,n);
    index(3) = 2*elements(2,n) - 1;
    index(4) = 2*elements(2,n);

    for ii = 1:4
        for jj = 1:4
            kk(index(ii), index(jj)) = kk(index(ii), index(jj)) + k(ii,jj,n);
        end
    end
end

jj_orig = kk; % Cópia antes das condições de contorno (para reações)

% CONDIÇÕES DE CONTORNO
for i = 1:size(apoio, 2)
    kk(apoio(i), :) = 0;
    kk(:, apoio(i)) = 0;
    kk(apoio(i), apoio(i)) = 1.0;
end

% SOLUÇÃO: K*u = F  →  deslocamentos nodais [mm]
u  = kk \ F;
FR = jj_orig * u; % Reações nos apoios [N]

% DESLOCAMENTOS DOS DOFs LIVRES [mm]
u_livre = u(dofs_livres);

% DEFORMAÇÕES E TENSÕES
for n = 1:n_el
    index(1) = 2*elements(1,n) - 1;
    index(2) = 2*elements(1,n);
    index(3) = 2*elements(2,n) - 1;
    index(4) = 2*elements(2,n);

    d(:,n) = [c(n)  s(n)  0     0    ;
              0     0     c(n)  s(n) ] * u(index);

    eps(n) = (1/L(n)) * (-d(1,n) + d(2,n));
end

Sigma = E * eps; % Tensão em cada barra [MPa]

% PESO TOTAL [kg]
weight = 0;
for i = 1:n_el
    weight = weight + area(i) * L(i) * dens; % área[mm²] × comp[mm] × dens[kg/mm³]
end

end
