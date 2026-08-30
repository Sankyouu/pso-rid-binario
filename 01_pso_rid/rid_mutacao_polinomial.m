function x_novo = rid_mutacao_polinomial(x, x_l, x_u, eta)
% RID_MUTACAO_POLINOMIAL  Mutacao polinomial de um numero real.
%
% BLOCO 1 (PSO-RID) — auxiliar isolado e testavel.
%
% REFERENCIA:
%   [DF2011] Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%            Eq. (9), Secao 5.
%   Originalmente de: Deb, K. "Multi-Objective Optimization using
%   Evolutionary Algorithms", John Wiley & Sons, 2001 (ref. [9] de [DF2011]).
%
% -------------------------------------------------------------------------
% Eq. (9) [DF2011]:
%
%   x  <-  x + (x^(u) - x^(l)) * d_q
%
%   d_q = [ 2r + (1-2r) * ( (x^(u) - x)/(x^(u) - x^(l)) )^(eta+1) ]^(1/(eta+1)) - 1
%         se r < 0.5
%
%   d_q = 1 - [ 2(1-r) + (2r-1) * ( (x - x^(l))/(x^(u) - x^(l)) )^(eta+1) ]^(1/(eta+1))
%         caso contrario
%
%   onde r e um numero aleatorio em ]0,1[ e eta > 0 o indice de distribuicao
%   polinomial.
%
% [DF2011] Sec. 5: "The distribution index eta in Eq. (9) is assigned a random
% value, in the range of [25,45], in different runs of a problem."
%
% [D5] A versao anterior deste solver usava perturbacao GAUSSIANA simples
%      (randn * 5% do dominio) na busca local, em vez da Eq. (9).
%
% -------------------------------------------------------------------------
% ENTRADAS
%   x   : valor real atual
%   x_l : limite inferior do dominio
%   x_u : limite superior do dominio
%   eta : indice de distribuicao polinomial (> 0)
%
% SAIDA
%   x_novo : valor mutado, garantidamente dentro de [x_l, x_u]
%
% PROPRIEDADE: quanto MAIOR eta, menor a perturbacao (distribuicao mais
% concentrada em torno de x). Essa propriedade e verificada em
% 06_testes/TestMutacaoPolinomial.m.

if x_u <= x_l
    % Dominio degenerado: nada a mutar
    x_novo = x;
    return;
end

assert(eta > 0, 'rid_mutacao_polinomial:etaInvalido', ...
    'O indice de distribuicao eta deve ser positivo. Recebido: %g.', eta);

r     = rand;
delta = x_u - x_l;

if r < 0.5
    base = 2*r + (1 - 2*r) * ((x_u - x) / delta)^(eta + 1);
    d_q  = base^(1/(eta + 1)) - 1;
else
    base = 2*(1 - r) + (2*r - 1) * ((x - x_l) / delta)^(eta + 1);
    d_q  = 1 - base^(1/(eta + 1));
end

x_novo = x + delta * d_q;

% Salvaguarda numerica: a Eq. (9) e construida para manter x em [x_l, x_u],
% mas arredondamentos podem levar a excursoes de ordem 1e-16.
x_novo = max(x_l, min(x_u, x_novo));

end
