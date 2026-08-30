% =========================================================================
% PSO - Particle Swarm Optimization (Otimização por Enxame de Partículas)
% =========================================================================
% Minimiza uma função objetivo FCN(x) definida em arquivo externo.
% O algoritmo simula um "enxame" de partículas que se movem pelo espaço
% de busca, atraídas pela melhor posição que cada uma já visitou (lbest)
% e pela melhor posição já encontrada por qualquer partícula (gbest).
% =========================================================================

clc; clear all; close all;

% -------------------------------------------------------------------------
% PARÂMETROS DO PROBLEMA
% -------------------------------------------------------------------------
xmin = 0.1*ones(1,10);   % Limite inferior do espaço de busca (uma entrada por variável)
xmax = 20*ones(1,10);     % Limite superior do espaço de busca
nrvar = length(xmin); % Número de variáveis do problema (dimensão do espaço)

% -------------------------------------------------------------------------
% PARÂMETROS DO PSO
% -------------------------------------------------------------------------
lambda1 = 2.02;  % Coeficiente cognitivo: peso da atração pela melhor posição
                 % individual (memória da própria partícula)
lambda2 = 2.52;  % Coeficiente social: peso da atração pela melhor posição
                 % global (influência do enxame)
omega   = 0.4;   % Inércia: quanto da velocidade anterior é mantida a cada passo
                 % (valores menores = convergência mais rápida, menos exploração)
pop     = 20;    % Tamanho da população (número de partículas)
tol     = 1E-10; % Tolerância para critério de convergência
itermax = 300;   % Número máximo de iterações
rng(2);          % Semente do gerador aleatório (garante reprodutibilidade)

% =========================================================================
% VISUALIZAÇÃO DA SUPERFÍCIE DA FUNÇÃO OBJETIVO (apenas para 2 variáveis)
% =========================================================================
% Cria uma grade de pontos para plotar a superfície 3D da função
%[X1,X2] = meshgrid(xmin(1):0.2:xmax(1), xmin(2):0.2:xmax(2));
%for i = 1:length(X1)
%    for j = 1:length(X2)
%        Z(i,j) = FCN([X1(i,j); X2(i,j)]); % Avalia FCN em cada ponto da grade
%    end
%end
%figure(1)
%mesh(X1,X2,Z); hold on;  % Plota a superfície 3D
%view(140,45)              % Define o ângulo de visão do gráfico

figure(1)
title('Convergência do PSO')
xlabel('Iteração'); ylabel('Peso (lb)')

% =========================================================================
% INICIALIZAÇÃO DA POPULAÇÃO
% =========================================================================
gbest = 1E30;           % Melhor valor global inicializado com valor muito alto
k = 1;                  % Contador de iterações começa em 1
v = zeros(pop, nrvar);  % Velocidades iniciais = zero para todas as partículas

for i = 1:pop
    % Gera posição inicial aleatória para a partícula i dentro dos limites
    for j = 1:nrvar
        x(i,j) = xmin(j) + (xmax(j) - xmin(j))*rand;
    end

    y = FCN(x(i,:));     % Avalia a função objetivo na posição inicial

    lbest(i)    = y;     % Melhor valor individual: começa com o valor inicial
    xlbest(i,:) = x(i,:); % Melhor posição individual: começa com a posição inicial

    % Atualiza o melhor global se esta partícula for melhor
    if (y < gbest)
        gbest(k) = y;
        xgbest   = x(i,:); % Posição correspondente ao melhor global
    end
end

% Plota a posição do melhor indivíduo na iteração inicial
figure(1)
plot3(xgbest(1), xgbest(2), gbest(k), 'ro', 'MarkerSize', 13, 'MarkerFaceColor', 'r'); hold on

% =========================================================================
% LOOP PRINCIPAL DO PSO
% =========================================================================
flag = 0;
k = 2; % Começa em 2 porque k=1 foi usado na inicialização

while(flag == 0)

    % Herda o melhor global da iteração anterior como ponto de partida
    gbest(k) = gbest(k-1);

    % ---------------------------------------------------------------------
    % ATUALIZAÇÃO DE VELOCIDADE E POSIÇÃO DE CADA PARTÍCULA
    % ---------------------------------------------------------------------
    for i = 1:pop
        for j = 1:nrvar
            r1 = rand; % Número aleatório para componente cognitiva
            r2 = rand; % Número aleatório para componente social

            % Equação de atualização da velocidade:
            % v_nova = inércia*v_antiga
            %        + cognitivo * r1 * (melhor_posição_individual - posição_atual)
            %        + social    * r2 * (melhor_posição_global      - posição_atual)
            vnew(i,j) = omega   * v(i,j) ...
                      + lambda1 * r1 * (xlbest(i,j) - x(i,j)) ...
                      + lambda2 * r2 * (xgbest(j)   - x(i,j));

            % Atualiza a posição somando a nova velocidade
            xnew(i,j) = x(i,j) + vnew(i,j);

            % Aplica limites: se sair do espaço de busca, projeta na fronteira
            if (xnew(i,j) < xmin(j))
                xnew(i,j) = xmin(j);
            elseif (xnew(i,j) > xmax(j))
                xnew(i,j) = xmax(j);
            end
        end

        % Avalia a função objetivo na nova posição
        ynew = FCN(xnew(i,:));

        % -----------------------------------------------------------------
        % FIGURA 4: Trajetória das 3 primeiras partículas no espaço de busca
        % (útil para visualizar como as partículas se movem)
        % -----------------------------------------------------------------
        figure(4)
        zz = 1; % Índice da primeira partícula a rastrear
        if (i == zz || i == zz+1 || i == zz+2)
            if (i == zz)
                z1(1:2, k-1) = [xnew(zz,1), xnew(zz,2)];
                plot(z1(1,:), z1(2,:), 'b'); hold on; % Partícula 1: azul
            elseif (i == zz+1)
                z2(1:2, k-1) = [xnew(zz+1,1), xnew(zz+1,2)];
                plot(z2(1,:), z2(2,:), 'r'); hold on; % Partícula 2: vermelho
            elseif (i == zz+2)
                z3(1:2, k-1) = [xnew(zz+2,1), xnew(zz+2,2)];
                plot(z3(1,:), z3(2,:), 'g'); hold on; % Partícula 3: verde
            end
            axis([xmin(1) xmax(1) xmin(2) xmax(2)])
            axis('equal')
        end

        % -----------------------------------------------------------------
        % ATUALIZAÇÃO DOS MELHORES (individual e global)
        % -----------------------------------------------------------------

        % Se a nova posição é melhor que o histórico individual, atualiza
        if (ynew < lbest(i))
            lbest(i)    = ynew;
            xlbest(i,:) = xnew(i,:);
        end

        % Se a nova posição é melhor que o melhor global atual, atualiza
        if (ynew < gbest(k))
            gbest(k) = ynew;
            xgbest   = xnew(i,:);
        end
    end

    % ---------------------------------------------------------------------
    % FIGURA 3: Convergência e localização do melhor resultado por iteração
    % ---------------------------------------------------------------------
    figure(3)
    subplot(2,1,1)
    plot(k, gbest(k), 'b*'); hold on;       % Evolução do melhor valor global
    axis([0 Inf 0 Inf])

    subplot(2,1,2)
    plot(xgbest(1), xgbest(2), 'b*'); hold on; % Posição do melhor global
    plot(1, 1, 'ro')                            % Marca o ótimo conhecido (1,1)
    axis([xmin(1) xmax(1) xmin(2) xmax(2)])

    % Se houve melhora, marca o novo melhor na superfície 3D (Figura 1)
    if (gbest(k) < gbest(k-1))
        figure(1)
        plot3(xgbest(1), xgbest(2), gbest(k), 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r'); hold on
    end

    % ---------------------------------------------------------------------
    % CRITÉRIOS DE PARADA
    % ---------------------------------------------------------------------

    % Critério 1: número máximo de iterações atingido
    if (k >= itermax)
        flag = 1;
    end

    % Critério 2: convergência por estagnação (diferença entre janelas de
    % iterações consecutivas abaixo da tolerância) — atualmente desativado
    if (k > 11)
        norm = sum(gbest(k-9:k-5)) - sum(gbest(k-4:k));
        if (norm < tol)
            % flag = 1;  % Descomente para ativar este critério
        end
    end

    % Avança para a próxima iteração e atualiza posições/velocidades
    k = k + 1;
    x = xnew;
    v = vnew;

end

% =========================================================================
% RESULTADOS FINAIS
% =========================================================================
disp('k = ')
disp(k-1)         % Número total de iterações executadas

disp('norm = ')
disp(norm)        % Último valor da métrica de convergência

disp('gbest = ')
disp(gbest(k-1))  % Melhor valor da função objetivo encontrado

disp('xgbest = ')
disp(xgbest)      % Posição (vetor de variáveis) do melhor resultado