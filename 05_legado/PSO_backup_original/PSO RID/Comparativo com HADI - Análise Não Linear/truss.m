function [weight, Sigma, u_livre] = truss(area)

% Análise GEOMETRICAMENTE NÃO LINEAR de treliça plana de 10 barras.
%
%
% Método: Newton-Raphson incremental (Eq. 13–18 do artigo)
%   [KT] = [KE] + [KG]                            (Eq. 13)
%   [KT]_i {Δu} = {ΔF}                            (Eq. 18)
%
% [KE] — matriz de rigidez elástica linear (idêntica à análise linear)
% [KG] — matriz de rigidez geométrica (Eq. 14), função das forças axiais
%         e ângulos de orientação atualizados a cada iteração
%
% Unidades:
%   Comprimento : mm
%   Força       : N
%   Tensão      : MPa (= N/mm²)
%   Massa       : kg
%
% USO:
%   [weight, Sigma, u_livre] = truss([19355,65,16129,7742,65,645,5161,12903,12903,65])
%
% REFERÊNCIA (Tabela 1, Case 2 discreto): peso = 2325.2 kg

%  PARÂMETROS DE CONTROLE DO NEWTON-RAPHSON INCREMENTAL
n_inc    = 10;       % Número de incrementos de carga (mais incrementos → mais estável)
max_iter = 50;       % Máximo de iterações N-R por incremento
tol      = 1e-6;     % Tolerância de convergência (norma relativa do resíduo)

%  PROPRIEDADES DO MATERIAL
E    = 6.89e4;       % Módulo de elasticidade [MPa]
dens = 2.770e-6;     % Densidade volumétrica  [kg/mm³]

%  GEOMETRIA: POSIÇÃO INICIAL DOS NÓS [mm]
%
%   nó 5 ------- nó 3 ------- nó 1
%        \       |       \       |
%         \      |        \      |
%   nó 6 ------- nó 4 ------- nó 2
%
%      ←── 9144 mm ──→←── 9144 mm ──→
%      ↑ 9144 mm ↓
%
%   nó    1  2  3  4  5  6
nodes0 = [2  2  1  1  0  0 ;   % coord x  [×9144 mm]
          1  0  1  0  1  0 ] * 9144;  % coord y

%  TOPOLOGIA: CONECTIVIDADE DAS BARRAS
%   el      1  2  3  4  5  6  7  8  9 10
elements = [5  3  6  4  3  1  5  6  3  4 ;   % nó inicial
            3  1  4  2  4  2  4  3  2  1];   % nó final

%  DIMENSÕES DO SISTEMA
n_el    = size(elements, 2);   % 10 barras
n_nodes = size(nodes0,   2);   % 6 nós
s_dof   = n_nodes * 2;         % 12 DOFs totais (2 por nó: x e y)

%  APOIOS: nós 5 e 6 completamente fixos
apoios      = [9 10 11 12];
todos_dofs  = 1:s_dof;
dofs_livres = setdiff(todos_dofs, apoios);

%  CARGA TOTAL APLICADA [N]
F_total = zeros(s_dof, 1);
F_total(4)  = -445400;   % nó 2, direção y
F_total(8)  = -445400;   % nó 4, direção y

%  INICIALIZAÇÃO DO NEWTON-RAPHSON INCREMENTAL
u = zeros(s_dof, 1);   % Deslocamentos acumulados (globais) [mm]
P_axial = zeros(n_el, 1);   % Forças axiais acumuladas nas barras [N]
nodes_cur = nodes0;          % Coordenadas atualizadas (cofiguração corrente)

% Comprimentos e ângulos iniciais (usados no cálculo do peso)
L0 = zeros(n_el, 1);
for n = 1:n_el
    ni = elements(1,n);  nj = elements(2,n);
    dx = nodes0(1,nj) - nodes0(1,ni);
    dy = nodes0(2,nj) - nodes0(2,ni);
    L0(n) = sqrt(dx^2 + dy^2);
end

%  LOOP DE INCREMENTOS DE CARGA
dF = F_total / n_inc;   % Incremento de carga em cada passo

for inc = 1:n_inc

    % --- Resíduo inicial do incremento: forças não equilibradas ---
    % No início de cada incremento o resíduo é simplesmente o incremento
    % de carga externo (as forças internas já equilibraram o passo anterior)
    R = dF;   % {R} = {ΔF} - ({F_int_novo} - {F_int_ant}) → inicia como ΔF

    % -----------------------------------------------------------------------
    %  ITERAÇÕES DE NEWTON-RAPHSON
    %  Eq. 18:  [KT]_i {δu} = {R}_i
    %  onde {R}_i = {F_ext} - {F_int} é o vetor de resíduo (forças não equil.)
    % -----------------------------------------------------------------------
    for iter = 1:max_iter

        % ----------------------------------------------------------------
        %  MONTAGEM DE [KT] = [KE] + [KG]  (Eq. 13 e 14)
        %  Usa a configuração corrente (nodes_cur) e as forças axiais (P_axial)
        % ----------------------------------------------------------------
        KT = zeros(s_dof, s_dof);

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            % Coordenadas correntes do elemento
            xi = nodes_cur(1,ni);  yi = nodes_cur(2,ni);
            xj = nodes_cur(1,nj);  yj = nodes_cur(2,nj);

            dx = xj - xi;
            dy = yj - yi;
            Ln = sqrt(dx^2 + dy^2);   % Comprimento corrente

            % Cossenos diretores correntes (ângulo de orientação atualizado)
            lam = dx / Ln;   % λ = cos(φ)  — Eq. 14
            mu  = dy / Ln;   % μ = sin(φ)  — Eq. 14

            % --- Matriz de rigidez elástica [KE] do elemento ---
            % Idêntica à análise linear, mas com geometria corrente
            ke = (area(n) * E / Ln) * ...
                [ lam*lam   lam*mu  -lam*lam  -lam*mu ;
                  lam*mu    mu*mu   -lam*mu   -mu*mu  ;
                 -lam*lam  -lam*mu   lam*lam   lam*mu ;
                 -lam*mu   -mu*mu    lam*mu    mu*mu  ];

            % --- Matriz de rigidez geométrica [KG] do elemento (Eq. 14) ---
            % Captura os efeitos de segunda ordem devidos à carga axial P
            % KG depende da força axial corrente P_axial(n) e do comprimento Ln
            Pn = P_axial(n);   % Força axial corrente nesta barra [N]

            kg = (Pn / Ln) * ...
                [  mu*mu   -lam*mu  -mu*mu    lam*mu ;
                  -lam*mu   lam*lam  lam*mu  -lam*lam;
                  -mu*mu    lam*mu   mu*mu   -lam*mu ;
                   lam*mu  -lam*lam -lam*mu   lam*lam];

            % --- Índices globais dos DOFs deste elemento ---
            idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];

            % --- Adição ao [KT] global ---
            KT(idx, idx) = KT(idx, idx) + ke + kg;
        end

        % ----------------------------------------------------------------
        %  CONDIÇÕES DE CONTORNO: zera linhas/colunas dos DOFs fixos
        %  e coloca 1 na diagonal → sistema compatível com restrições
        % ----------------------------------------------------------------
        KT_bc = KT;
        R_bc  = R;
        for dof = apoios
            KT_bc(dof, :) = 0;
            KT_bc(:, dof) = 0;
            KT_bc(dof, dof) = 1.0;
            R_bc(dof) = 0.0;
        end

        % ----------------------------------------------------------------
        %  SOLUÇÃO: δu = [KT]^{-1} {R}
        %  Deslocamentos incrementais desta iteração N-R
        % ----------------------------------------------------------------
        du = KT_bc \ R_bc;

        % ----------------------------------------------------------------
        %  ATUALIZAÇÃO DOS DESLOCAMENTOS E CONFIGURAÇÃO
        % ----------------------------------------------------------------
        u = u + du;

        % Atualiza as coordenadas correntes dos nós
        for nd = 1:n_nodes
            nodes_cur(1,nd) = nodes0(1,nd) + u(2*nd-1);
            nodes_cur(2,nd) = nodes0(2,nd) + u(2*nd);
        end

        % ----------------------------------------------------------------
        %  ATUALIZAÇÃO DAS FORÇAS AXIAIS (configuração corrente)
        %  Baseado na deformação axial de cada barra com a geometria nova
        % ----------------------------------------------------------------
        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            % Comprimento inicial da barra
            dx0 = nodes0(1,nj) - nodes0(1,ni);
            dy0 = nodes0(2,nj) - nodes0(2,ni);
            L_ini = sqrt(dx0^2 + dy0^2);

            % Comprimento corrente da barra
            dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
            dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
            L_cor = sqrt(dxc^2 + dyc^2);

            % Deformação de engenharia rotacionada (Crisfield Eq. 3.103/3.104)
            % ε = (L_cor - L_ini) / L_ini
            eps_n = (L_cor - L_ini) / L_ini;

            % Força axial = EA × ε
            P_axial(n) = area(n) * E * eps_n;
        end

        % ----------------------------------------------------------------
        %  CÁLCULO DO RESÍDUO: R = F_ext_acumulado - F_int_corrente
        %  F_ext_acumulado = inc × dF   (carga aplicada até este incremento)
        %  F_int é calculado das forças axiais via equilíbrio nodal
        % ----------------------------------------------------------------
        F_ext_acc = inc * dF;   % Carga total aplicada até agora
        F_int     = zeros(s_dof, 1);

        for n = 1:n_el
            ni = elements(1,n);
            nj = elements(2,n);

            dxc = nodes_cur(1,nj) - nodes_cur(1,ni);
            dyc = nodes_cur(2,nj) - nodes_cur(2,ni);
            Ln  = sqrt(dxc^2 + dyc^2);

            lam = dxc / Ln;
            mu  = dyc / Ln;

            % Contribuição desta barra ao vetor de forças internas
            % q_int = P × [−λ, −μ, +λ, +μ]^T  (equilíbrio do elemento)
            idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
            felem = P_axial(n) * [-lam; -mu; lam; mu];
            F_int(idx) = F_int(idx) + felem;
        end

        % Resíduo = forças externas - forças internas (desequilíbrio)
        R = F_ext_acc - F_int;
        R(apoios) = 0;   % Não há resíduo nos DOFs restritos

        % ----------------------------------------------------------------
        %  CRITÉRIO DE CONVERGÊNCIA
        %  Norma relativa do resíduo em relação à carga aplicada
        % ----------------------------------------------------------------
        norma_R  = norm(R);
        norma_F  = norm(F_ext_acc);
        if norma_F < 1e-12, norma_F = 1.0; end   % Evita divisão por zero

        if norma_R / norma_F < tol
            break   % Convergiu: passa para o próximo incremento
        end

    end  % fim das iterações N-R
end  % fim dos incrementos

%  PÓS-PROCESSAMENTO

% --- Deslocamentos dos DOFs livres [mm] ---
u_livre = u(dofs_livres);

% --- Tensões finais nas barras [MPa] ---
Sigma = P_axial ./ area(:);   % σ = P / A  (tração > 0, compressão < 0)

% --- Peso total [kg] ---
% Usa comprimentos iniciais (L0), pois o peso é propriedade geométrica do projeto
weight = sum(area(:) .* L0 .* dens);

end