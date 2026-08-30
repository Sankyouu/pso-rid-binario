function [v_next, x_next] = rid_velocidade_binaria(x_atual, tendencia, pm)
% RID_VELOCIDADE_BINARIA  Mapeia a velocidade real para {-1,0,+1} numa dimensao binaria.
%
% BLOCO 1 (PSO-RID) — auxiliar isolado e testavel.
%
% REFERENCIA:
%   [DF2011] Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%            Secao 4.3, Eq. (7), Eq. (8), Tabelas 1 e 2.
%
% -------------------------------------------------------------------------
% Eq. (7) [DF2011]: transicoes de posicao binaria admissiveis
%
%   x^(i,t+1) = 0  se (x^(i,t), v^(i,t+1)) = (0,0) ou (1,-1)
%   x^(i,t+1) = 1  se (x^(i,t), v^(i,t+1)) = (0,1) ou (1,0)
%
%   Ou seja:   x = 0  =>  v pertence a { 0, +1}     (Tabela 1)
%              x = 1  =>  v pertence a { 0, -1}     (Tabela 2)
%
% [D3] Esta assimetria e o ponto central: a versao anterior deste solver
%      usava a regra SIMETRICA v em {-1,0,+1} independentemente do bit atual,
%      saturando a posicao depois. Isso corrompia a memoria de inercia
%      (w * v no passo seguinte) em varios casos das Tabelas 1 e 2.
%
% -------------------------------------------------------------------------
% Eq. (8) [DF2011]: mapeamento do valor real de velocidade para o discreto
%
%   v > 0                          =>  v = +1
%   v < 0                          =>  v = -1
%   v >= 0, ou sem relacao na Tab.1 =>  v = 0 ou +1
%   v <= 0, ou sem relacao na Tab.2 =>  v = 0 ou -1
%
% A ambiguidade ("0 ou 1" / "0 ou -1") e resolvida pela probabilidade de
% mutacao p_m (coluna 7 das Tabelas 1 e 2): com probabilidade p_m adota-se o
% valor ALTERNATIVO ao mapeado deterministicamente; caso contrario, o mapeado.
%
% [D4] A versao anterior sorteava uniformemente em {-1,0,+1} sem condicionar
%      ao bit atual, desperdicando parte dos sorteios na saturacao posterior.
%
% [DF2011] Sec. 4.3 sobre p_m: "Extensive empirical studies have shown that
% better results are obtained for p_m <= 15%."
%
% -------------------------------------------------------------------------
% ENTRADAS
%   x_atual   : bit atual da dimensao (0 ou 1)
%   tendencia : valor REAL da velocidade calculado pela Eq. (1)
%               v = w*v + c1*r1*(pbest - x) + c2*r2*(gbest - x)
%   pm        : probabilidade de mutacao (adocao do valor alternativo)
%
% SAIDAS
%   v_next : velocidade discreta resultante (-1, 0 ou +1)
%   x_next : nova posicao binaria, x_atual + v_next (garantidamente 0 ou 1)

validateattributes(x_atual, {'numeric'}, {'scalar'});
assert(x_atual == 0 || x_atual == 1, ...
    'rid_velocidade_binaria:posicaoInvalida', ...
    'x_atual deve ser 0 ou 1. Recebido: %g.', x_atual);

if x_atual == 0
    % ---- Tabela 1 [DF2011]: valores admissiveis de v sao 0 ou +1 ----
    if tendencia > 0
        v_det = 1;
    elseif tendencia < 0
        v_det = 0;    % -1 e inadmissivel para x=0 (Eq. 7)
    else
        v_det = double(rand < 0.5);   % empate exato: 50% cada
    end
    v_alt = 1 - v_det;                % alternativa dentro de {0, +1}

else
    % ---- Tabela 2 [DF2011]: valores admissiveis de v sao 0 ou -1 ----
    if tendencia < 0
        v_det = -1;
    elseif tendencia > 0
        v_det = 0;    % +1 e inadmissivel para x=1 (Eq. 7)
    else
        v_det = -double(rand < 0.5);  % empate exato: 50% cada
    end
    v_alt = -1 - v_det;               % alternativa dentro de {0, -1}
end

% Randomizacao pela probabilidade de mutacao p_m
% (coluna 7 das Tabelas 1 e 2 [DF2011])
if rand < pm
    v_next = v_alt;
else
    v_next = v_det;
end

% Eq. (2)/(7) [DF2011]: por construcao das transicoes admissiveis, a soma
% resulta sempre em {0,1} — nao ha necessidade de saturacao.
x_next = x_atual + v_next;

end
