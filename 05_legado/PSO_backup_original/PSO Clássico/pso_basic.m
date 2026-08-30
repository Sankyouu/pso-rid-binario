clc; clear all; close all;

% PARÂMETROS DO PROBLEMA
xmin = 0.1*ones(1,10);
xmax = 20*ones(1,10);
nrvar = length(xmin);

% PARÂMETROS DO PSO
lambda1 = 2.02;
lambda2 = 2.52;
omega   = 0.4;
pop     = 20;
tol     = 1E-10;
itermax = 300;
rng(2);

% VISUALIZAÇÃO DA SUPERFÍCIE DA FUNÇÃO OBJETIVO (apenas para 2 variáveis)
%[X1,X2] = meshgrid(xmin(1):0.2:xmax(1), xmin(2):0.2:xmax(2));
%for i = 1:length(X1)
%    for j = 1:length(X2)
%        Z(i,j) = FCN([X1(i,j); X2(i,j)]);
%    end
%end
%figure(1)
%mesh(X1,X2,Z); hold on;
%view(140,45)

figure(1)
title('Convergência do PSO')
xlabel('Iteração'); ylabel('Peso (lb)')

% INICIALIZAÇÃO DA POPULAÇÃO
gbest = 1E30;
k = 1;
v = zeros(pop, nrvar);

for i = 1:pop
    for j = 1:nrvar
        x(i,j) = xmin(j) + (xmax(j) - xmin(j))*rand;
    end

    y = FCN(x(i,:));

    lbest(i)    = y;
    xlbest(i,:) = x(i,:);

    if (y < gbest)
        gbest(k) = y;
        xgbest   = x(i,:);
    end
end

figure(1)
plot3(xgbest(1), xgbest(2), gbest(k), 'ro', 'MarkerSize', 13, 'MarkerFaceColor', 'r'); hold on

% LOOP PRINCIPAL DO PSO
flag = 0;
k = 2; % Começa em 2 porque k=1 foi usado na inicialização

while(flag == 0)

    gbest(k) = gbest(k-1);

    % ATUALIZAÇÃO DE VELOCIDADE E POSIÇÃO DE CADA PARTÍCULA
    for i = 1:pop
        for j = 1:nrvar
            r1 = rand;
            r2 = rand;

            vnew(i,j) = omega   * v(i,j) ...
                      + lambda1 * r1 * (xlbest(i,j) - x(i,j)) ...
                      + lambda2 * r2 * (xgbest(j)   - x(i,j));

            xnew(i,j) = x(i,j) + vnew(i,j);

            if (xnew(i,j) < xmin(j))
                xnew(i,j) = xmin(j);
            elseif (xnew(i,j) > xmax(j))
                xnew(i,j) = xmax(j);
            end
        end

        ynew = FCN(xnew(i,:));

        % FIGURA 4: Trajetória das 3 primeiras partículas no espaço de busca
        figure(4)
        zz = 1;
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

        % ATUALIZAÇÃO DOS MELHORES (individual e global)

        if (ynew < lbest(i))
            lbest(i)    = ynew;
            xlbest(i,:) = xnew(i,:);
        end

        if (ynew < gbest(k))
            gbest(k) = ynew;
            xgbest   = xnew(i,:);
        end
    end

    % FIGURA 3: Convergência e localização do melhor resultado por iteração
    figure(3)
    subplot(2,1,1)
    plot(k, gbest(k), 'b*'); hold on;
    axis([0 Inf 0 Inf])

    subplot(2,1,2)
    plot(xgbest(1), xgbest(2), 'b*'); hold on;
    plot(1, 1, 'ro')
    axis([xmin(1) xmax(1) xmin(2) xmax(2)])

    if (gbest(k) < gbest(k-1))
        figure(1)
        plot3(xgbest(1), xgbest(2), gbest(k), 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r'); hold on
    end

    % CRITÉRIOS DE PARADA

    % Critério 1: número máximo de iterações atingido
    if (k >= itermax)
        flag = 1;
    end

    % Critério 2: convergência por estagnação (diferença entre janelas de
    % iterações consecutivas abaixo da tolerância)
    if (k > 11)
        norm = sum(gbest(k-9:k-5)) - sum(gbest(k-4:k));
        if (norm < tol)
             flag = 1;
        end
    end

    % Avança para a próxima iteração e atualiza posições/velocidades
    k = k + 1;
    x = xnew;
    v = vnew;

end

% RESULTADOS FINAIS
disp('k = ')
disp(k-1)         % Número total de iterações executadas

disp('norm = ')
disp(norm)        % Último valor da métrica de convergência

disp('gbest = ')
disp(gbest(k-1))  % Melhor valor da função objetivo encontrado

disp('xgbest = ')
disp(xgbest)      % Posição (vetor de variáveis) do melhor resultado

[weight_final, Sigma_final] = truss(xgbest);
disp('Peso real (sem penalidade) = ')
disp(weight_final)