function y = FCN(area)

    % PARÂMETRO DE RESTRIÇÃO
    sigma_max = 25;
    [weight, Sigma] = truss(area);

    % FUNÇÃO DE PENALIZAÇÃO POR VIOLAÇÃO DE TENSÃO
    penalty = 0;
    for i = 1:length(Sigma)
        if abs(Sigma(i)) > sigma_max
            penalty = penalty + 1e6 * (abs(Sigma(i)) - sigma_max)^2;
        end
    end

    % VALOR FINAL DA FUNÇÃO OBJETIVO
    y = weight + penalty;

end