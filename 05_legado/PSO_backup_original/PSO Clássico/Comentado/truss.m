% =========================================================================
% TRUSS.M - Análise Estrutural de Treliça por Elementos Finitos
% =========================================================================
% Resolve uma treliça plana de 10 barras e 6 nós usando o Método dos
% Elementos Finitos (MEF). Calcula o peso total e as tensões em cada barra
% para um dado vetor de áreas transversais.
%
% Entrada:
%   area   - vetor (1x10) com a área transversal de cada barra [in²]
%
% Saídas:
%   weight - peso total da estrutura [lb]
%   Sigma  - vetor (1x10) com a tensão normal em cada barra [ksi]
%            positivo = tração, negativo = compressão
% =========================================================================

function [weight, Sigma] = truss(area)

% -------------------------------------------------------------------------
% PROPRIEDADES DO MATERIAL
% -------------------------------------------------------------------------
E    = 29000; % Módulo de elasticidade (Young) [ksi] — equivale a ~200 GPa (aço)
dens = 0.1;   % Densidade do material [lb/in³]

% -------------------------------------------------------------------------
% GEOMETRIA: POSIÇÃO DOS NÓS
% -------------------------------------------------------------------------
% Cada coluna representa um nó: [x; y] em polegadas.
% O fator 360 converte as coordenadas normalizadas para polegadas.
% Layout (numeração dos nós):
%
%   5 ------- 3 ------- 1
%   |    \    |    \    |
%   6 ------- 4 ------- 2
%
%       360in     360in
%
nodes = 360*[2 2 1 1 0 0 ;   % coordenadas x dos nós 1..6
             1 0 1 0 1 0];   % coordenadas y dos nós 1..6

% -------------------------------------------------------------------------
% TOPOLOGIA: CONECTIVIDADE DOS ELEMENTOS
% -------------------------------------------------------------------------
% Cada coluna define um elemento (barra): [nó_i; nó_j]
% O elemento conecta o nó_i ao nó_j.
%
%   El:     1  2  3  4  5  6  7  8  9 10
elements = [5  3  6  4  3  1  5  6  3  4 ;   % nó inicial de cada barra
            3  1  4  2  4  2  4  3  2  1];   % nó final   de cada barra

% -------------------------------------------------------------------------
% DIMENSÕES DO SISTEMA
% -------------------------------------------------------------------------
n_el    = size(elements, 2); % Número de elementos (barras) = 10
ndof    = 2;                 % Graus de liberdade por nó (x e y)
n_nodes = size(nodes, 2);    % Número de nós = 6
s_dof   = n_nodes * ndof;    % Total de DOFs do sistema = 12

% -------------------------------------------------------------------------
% INICIALIZAÇÃO DA MATRIZ DE RIGIDEZ GLOBAL
% -------------------------------------------------------------------------
kk = zeros(s_dof, s_dof); % Matriz de rigidez global (12x12), zerada

% -------------------------------------------------------------------------
% CONDIÇÕES DE CONTORNO: APOIOS
% -------------------------------------------------------------------------
% DOFs restringidos pelos apoios (nós fixos ou roletes).
% Índices referem-se à numeração global dos graus de liberdade:
%   DOF 2k-1 = deslocamento horizontal do nó k
%   DOF 2k   = deslocamento vertical   do nó k
% Nós 5 e 6 (DOFs 9,10,11,12) estão fixos (apoio duplo).
apoio = [9 10 11 12];

% -------------------------------------------------------------------------
% VETOR DE FORÇAS EXTERNAS
% -------------------------------------------------------------------------
% Forças nodais no sistema global [kips], ordenadas por DOF.
% DOF 4  (vertical do nó 2): -100 kips (força para baixo)
% DOF 8  (vertical do nó 4): -100 kips (força para baixo)
% Demais DOFs: sem força aplicada
F = [0; 0; 0; -100; 0; 0; 0; -100; 0; 0; 0; 0];

% =========================================================================
% MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL
% =========================================================================
for n = 1:n_el
    i = elements(1,n); % Nó inicial do elemento n
    j = elements(2,n); % Nó final   do elemento n

    % Comprimento da barra (distância euclidiana entre os nós)
    L(n) = sqrt((nodes(1,j)-nodes(1,i))^2 + (nodes(2,j)-nodes(2,i))^2);

    % Cossenos diretores: projeção da barra nos eixos x e y
    c(n) = (nodes(1,j) - nodes(1,i)) / L(n); % cos(theta)
    s(n) = (nodes(2,j) - nodes(2,i)) / L(n); % sin(theta)

    % -----------------------------------------------------------------
    % MATRIZ DE RIGIDEZ LOCAL DO ELEMENTO (coordenadas globais)
    % -----------------------------------------------------------------
    % Para uma barra de treliça, a rigidez local 4x4 em coord. globais é:
    %   k = (A*E/L) * [T]^T * [1 -1; -1 1] * [T]
    % onde [T] é a matriz de transformação de coordenadas.
    % O resultado expandido é a matriz abaixo:
    k(:,:,n) = (area(n)*E/L(n)) * ...
        [ c(n)*c(n)  c(n)*s(n) -c(n)*c(n) -c(n)*s(n) ;
          c(n)*s(n)  s(n)*s(n) -c(n)*s(n) -s(n)*s(n) ;
         -c(n)*c(n) -c(n)*s(n)  c(n)*c(n)  c(n)*s(n) ;
         -c(n)*s(n) -s(n)*s(n)  c(n)*s(n)  s(n)*s(n) ];

    % -----------------------------------------------------------------
    % MAPEAMENTO: DOFs locais → DOFs globais
    % -----------------------------------------------------------------
    % Cada elemento tem 4 DOFs locais (2 por nó).
    % O índice global do DOF horizontal do nó k é 2k-1, e o vertical é 2k.
    index(1) = 2*elements(1,n) - 1; % DOF x do nó i
    index(2) = 2*elements(1,n);     % DOF y do nó i
    index(3) = 2*elements(2,n) - 1; % DOF x do nó j
    index(4) = 2*elements(2,n);     % DOF y do nó j

    % -----------------------------------------------------------------
    % ASSEMBLÉIA: soma a rigidez local na posição correta da matriz global
    % -----------------------------------------------------------------
    for i = 1:4
        ii = index(i);
        for j = 1:4
            jj = index(j);
            kk(ii,jj) = kk(ii,jj) + k(i,j,n);
        end
    end
end

% Guarda uma cópia da matriz de rigidez antes de aplicar as CCs
% (necessária depois para calcular as reações de apoio)
jj = kk;

% =========================================================================
% APLICAÇÃO DAS CONDIÇÕES DE CONTORNO (APOIOS)
% =========================================================================
% Método da linha/coluna zerada com 1 na diagonal:
% Para cada DOF restringido, zeramos sua linha e coluna na matriz de rigidez
% e colocamos 1 na diagonal. Isso força o deslocamento correspondente a ser
% zero na solução do sistema linear.
for i = 1:size(apoio, 2)
    kk(apoio(i), :) = 0;          % Zera a linha do DOF restringido
    kk(:, apoio(i)) = 0;          % Zera a coluna do DOF restringido
    kk(apoio(i), apoio(i)) = 1.;  % Coloca 1 na diagonal
end

% =========================================================================
% SOLUÇÃO DO SISTEMA LINEAR: K * u = F
% =========================================================================
u = kk \ F; % Vetor de deslocamentos nodais globais [in]
            % u(2k-1) = deslocamento horizontal do nó k
            % u(2k)   = deslocamento vertical   do nó k

% Calcula as forças de reação nos apoios usando a matriz de rigidez original
FR = jj * u; % FR nos DOFs livres será ~0; nos apoios dá as reações [kips]

% =========================================================================
% CÁLCULO DAS DEFORMAÇÕES E TENSÕES EM CADA BARRA
% =========================================================================
for n = 1:n_el
    i = elements(1,n);
    j = elements(2,n);

    % Mapeamento de DOFs locais para globais (mesmo que antes)
    index(1) = 2*elements(1,n) - 1;
    index(2) = 2*elements(1,n);
    index(3) = 2*elements(2,n) - 1;
    index(4) = 2*elements(2,n);

    % -----------------------------------------------------------------
    % TRANSFORMAÇÃO PARA COORDENADAS LOCAIS
    % -----------------------------------------------------------------
    % Projeta os deslocamentos globais dos nós do elemento no eixo da barra.
    % d(1,n) = deslocamento axial do nó i (projeção de u_i na direção da barra)
    % d(2,n) = deslocamento axial do nó j (projeção de u_j na direção da barra)
    d(:,n) = [c(n)  s(n)  0     0    ;
              0     0     c(n)  s(n) ] * u(index);

    % -----------------------------------------------------------------
    % DEFORMAÇÃO AXIAL (strain) da barra n
    % -----------------------------------------------------------------
    % eps = (d_j - d_i) / L  →  alongamento relativo
    % d(1,n) é o deslocamento do nó inicial, d(2,n) do nó final
    eps(n) = (1/L(n)) * (-d(1,n) + d(2,n));
end

% =========================================================================
% TENSÕES E PESO
% =========================================================================

% Tensão normal em cada barra pela Lei de Hooke: sigma = E * eps [ksi]
% Positivo = tração, negativo = compressão
Sigma = E * eps;

% Peso total da estrutura: soma de (área × comprimento × densidade) por barra
weight = 0;
for i = 1:n_el
    Fbars(i) = Sigma(i) * area(i); % Força axial na barra i [kips]
    weight   = weight + area(i) * L(i) * dens; % Contribuição ao peso [lb]
end

% =========================================================================
% PLOTAGEM DA DEFORMADA (desativada por padrão)
% =========================================================================
% Descomente o bloco abaixo para visualizar a estrutura original (tracejado)
% e a deformada amplificada (linha sólida).
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