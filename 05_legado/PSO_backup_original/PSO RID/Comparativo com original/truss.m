function [weight, Sigma] = truss(area)

% PROPRIEDADES DO MATERIAL
E    = 29000; % Módulo de elasticidade (Young) [ksi]
dens = 0.1;   % Densidade do material [lb/in³]

% GEOMETRIA: POSIÇÃO DOS NÓS
%
%   5 ------- 3 ------- 1
%   |    \    |    \    |
%   6 ------- 4 ------- 2
%
%       360in     360in
%
nodes = 360*[2 2 1 1 0 0 ;   % coordenadas x dos nós 1..6
             1 0 1 0 1 0];   % coordenadas y dos nós 1..6

% TOPOLOGIA: CONECTIVIDADE DOS ELEMENTOS
%   El:     1  2  3  4  5  6  7  8  9 10
elements = [5  3  6  4  3  1  5  6  3  4 ;   % nó inicial de cada barra
            3  1  4  2  4  2  4  3  2  1];   % nó final   de cada barra

% DIMENSÕES DO SISTEMA
n_el    = size(elements, 2); % Número de elementos (barras) = 10
ndof    = 2;                 % Graus de liberdade por nó (x e y)
n_nodes = size(nodes, 2);    % Número de nós = 6
s_dof   = n_nodes * ndof;    % Total de DOFs do sistema = 12

% INICIALIZAÇÃO DA MATRIZ DE RIGIDEZ GLOBAL
kk = zeros(s_dof, s_dof); % Matriz de rigidez global (12x12), zerada

% CONDIÇÕES DE CONTORNO: APOIOS
apoio = [9 10 11 12];

% VETOR DE FORÇAS EXTERNAS
F = [0; 0; 0; -100; 0; 0; 0; -100; 0; 0; 0; 0];

% MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL
for n = 1:n_el
    i = elements(1,n); % Nó inicial do elemento n
    j = elements(2,n); % Nó final   do elemento n

    % Comprimento da barra (distância euclidiana entre os nós)
    L(n) = sqrt((nodes(1,j)-nodes(1,i))^2 + (nodes(2,j)-nodes(2,i))^2);

    % Cossenos diretores: projeção da barra nos eixos x e y
    c(n) = (nodes(1,j) - nodes(1,i)) / L(n); % cos(theta)
    s(n) = (nodes(2,j) - nodes(2,i)) / L(n); % sin(theta)

    % MATRIZ DE RIGIDEZ LOCAL DO ELEMENTO (coordenadas globais)
    k(:,:,n) = (area(n)*E/L(n)) * ...
        [ c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
          c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
         -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
         -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    % MAPEAMENTO: DOFs locais → DOFs globais
    index(1) = 2*elements(1,n) - 1; % DOF x do nó i
    index(2) = 2*elements(1,n);     % DOF y do nó i
    index(3) = 2*elements(2,n) - 1; % DOF x do nó j
    index(4) = 2*elements(2,n);     % DOF y do nó j

    % ASSEMBLÉIA: soma a rigidez local na posição correta da matriz global
    for i = 1:4
        ii = index(i);
        for j = 1:4
            jj = index(j);
            kk(ii,jj) = kk(ii,jj) + k(i,j,n);
        end
    end
end

jj = kk;

% APLICAÇÃO DAS CONDIÇÕES DE CONTORNO (APOIOS)
for i = 1:size(apoio, 2)
    kk(apoio(i), :) = 0;
    kk(:, apoio(i)) = 0;
    kk(apoio(i), apoio(i)) = 1.;
end

% SOLUÇÃO DO SISTEMA LINEAR: K * u = F
u = kk \ F;

FR = jj * u;

% CÁLCULO DAS DEFORMAÇÕES E TENSÕES EM CADA BARRA
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);

    % Mapeamento de DOFs locais para globais (mesmo que antes)
    index(1) = 2*elements(1,n) - 1;
    index(2) = 2*elements(1,n);
    index(3) = 2*elements(2,n) - 1;
    index(4) = 2*elements(2,n);

    % TRANSFORMAÇÃO PARA COORDENADAS LOCAIS
    d(:,n) = [c(n)  s(n)  0     0    ;
              0     0     c(n)  s(n) ] * u(index);

    % DEFORMAÇÃO AXIAL (strain) da barra n
    eps(n) = (1/L(n)) * (-d(1,n) + d(2,n));
end

% TENSÕES E PESO
Sigma = E * eps;

% Peso total da estrutura: soma de (área × comprimento × densidade) por barra
weight = 0;
for i = 1:n_el
    Fbars(i) = Sigma(i) * area(i);
    weight   = weight + area(i) * L(i) * dens;
end

end