% =========================================================================
% FCN.m - Função Objetivo para o PSO
% =========================================================================
% Interface entre o PSO e o problema estrutural da treliça.
% Recebe um vetor de áreas das barras, calcula o peso total da estrutura
% e penaliza soluções que violem o limite de tensão admissível.
%
% Entrada:
%   area  - vetor (1x10) com as áreas transversais de cada barra [in²]
%
% Saída:
%   y     - escalar com o valor da função objetivo (peso + penalidades) [lb]
% =========================================================================

function y = FCN(area)

    % ---------------------------------------------------------------------
    % PARÂMETRO DE RESTRIÇÃO
    % ---------------------------------------------------------------------
    sigma_max = 25; % Tensão admissível máxima [ksi]
                    % Nenhuma barra pode ter |tensão| maior que este valor.
                    % Ajuste conforme a norma ou material utilizado.

    % ---------------------------------------------------------------------
    % CHAMADA DO SOLVER ESTRUTURAL
    % ---------------------------------------------------------------------
    % Chama a função truss.m que monta e resolve o sistema de elementos
    % finitos da treliça, retornando:
    %   weight - peso total da estrutura [lb]
    %   Sigma  - vetor (1x10) com a tensão em cada barra [ksi]
    [weight, Sigma] = truss(area);

    % ---------------------------------------------------------------------
    % FUNÇÃO DE PENALIZAÇÃO POR VIOLAÇÃO DE TENSÃO
    % ---------------------------------------------------------------------
    % O PSO é um método irrestrito (não trata restrições diretamente).
    % Para forçar o respeito ao limite de tensão, adicionamos uma penalidade
    % quadrática ao peso: quanto maior a violação, maior o valor de y,
    % fazendo o PSO rejeitar naturalmente essas soluções.
    penalty = 0;
    for i = 1:length(Sigma)
        if abs(Sigma(i)) > sigma_max
            % Penalidade proporcional ao quadrado da violação.
            % O fator 1e6 garante que qualquer violação domine o peso,
            % tornando a solução inviável muito pior que qualquer viável.
            penalty = penalty + 1e6 * (abs(Sigma(i)) - sigma_max)^2;
        end
    end

    % ---------------------------------------------------------------------
    % VALOR FINAL DA FUNÇÃO OBJETIVO
    % ---------------------------------------------------------------------
    % O PSO minimiza y. Para soluções viáveis (sem violação), y = weight.
    % Para soluções inviáveis, a penalidade empurra y para valores altos,
    % guiando o enxame de volta para a região viável do espaço de busca.
    y = weight + penalty;

end