function [weight, Sigma,u] = truss2(area)

% PROPRIEDADES DO MATERIAL
E    = 6.86e10;
dens = 2770;

% GEOMETRIA: POSIÇÃO DOS NÓS
nodes = 9144e-3*[2 2 1 1 0 0 ;
             1 0 1 0 1 0];

% TOPOLOGIA: CONECTIVIDADE DOS ELEMENTOS
%   El:     1  2  3  4  5  6  7  8  9 10
elements = [5  3  6  4  3  1  5  6  3  4 ;
            3  1  4  2  4  2  4  3  2  1];

% DIMENSÕES DO SISTEMA
n_el    = size(elements, 2);
ndof    = 2;
n_nodes = size(nodes, 2);
s_dof   = n_nodes * ndof;

% INICIALIZAÇÃO DA MATRIZ DE RIGIDEZ GLOBAL
kk = zeros(s_dof, s_dof);

% CONDIÇÕES DE CONTORNO: APOIOS
apoio = [9 10 11 12];

% VETOR DE FORÇAS EXTERNAS
F = [0; 0; 0; -445.4e3; 0; 0; 0; -445.4e3; 0; 0; 0; 0];

% MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);

    % Comprimento da barra (distância euclidiana entre os nós)
    L(n) = sqrt((nodes(1,j)-nodes(1,i))^2 + (nodes(2,j)-nodes(2,i))^2);

    % Cossenos diretores: projeção da barra nos eixos x e y
    c(n) = (nodes(1,j) - nodes(1,i)) / L(n);
    s(n) = (nodes(2,j) - nodes(2,i)) / L(n);

    % MATRIZ DE RIGIDEZ LOCAL DO ELEMENTO (coordenadas globais)
    k(:,:,n) = (area(n)*E/L(n)) * ...
        [ c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
          c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
         -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
         -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    % MAPEAMENTO: DOFs locais → DOFs globais
    index(1) = 2*elements(1,n) - 1;
    index(2) = 2*elements(1,n);
    index(3) = 2*elements(2,n) - 1;
    index(4) = 2*elements(2,n);

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

% PLOTAGEM DA DEFORMADA (desativada por padrão)
%
% amp = 10; % Fator de amplificação dos deslocamentos para visualização
%
% for i = 1:n_nodes
%     new_nodes(1,i) = nodes(1,i) + amp*u(2*i-1);
%     new_nodes(2,i) = nodes(2,i) + amp*u(2*i);
% end
%
% for i = 1:n_el
%     figure(1)
%     axis equal
%     axis(360*[-0.3 2.3 -0.3 1.3])
%     % Estrutura original (tracejado vermelho)
%     plot([nodes(1,elements(1,i))     nodes(1,elements(2,i))], ...
%          [nodes(2,elements(1,i))     nodes(2,elements(2,i))], '--rs')
%     hold on
%     % Estrutura deformada (linha azul)
%     plot([new_nodes(1,elements(1,i)) new_nodes(1,elements(2,i))], ...
%          [new_nodes(2,elements(1,i)) new_nodes(2,elements(2,i))], '-bo')
% end

end