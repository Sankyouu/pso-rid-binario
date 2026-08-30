function [vars, viol_estrutural] = rid_decodificar(vetor_hibrido, mapa, config, modo_discreto)
% RID_DECODIFICAR  Converte a particula hibrida (real + binaria) em valores de projeto.
%
% BLOCO 1 (PSO-RID) — auxiliar isolado e testavel.
%
% REFERENCIA:
%   [DF2011]  Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%
% -------------------------------------------------------------------------
% Eq. (6) [DF2011]:   x = sum_{i=1}^{B} 2^(B-i) * b_i
%   x   = valor inteiro da variavel
%   B   = numero de bits
%   b_i = i-esimo bit
% (codificacao binaria natural, bit mais significativo a esquerda)
% -------------------------------------------------------------------------
%
% [D1] VARIAVEL DISCRETA
%   [DF2011] Sec. 4.1 trata a discreta como inteira cujo valor representa o
%   INDICE do valor discreto real, com limite inferior 1 e limite superior N.
%   Portanto, no modo fiel ao artigo:
%         indice = x + 1          (x da Eq. 6, comecando em 0)
%
%   modo_discreto = 'datta'  (padrao, fiel ao artigo)
%       indice = x + 1, valido em 1..N. Indices > N sao ESTRUTURALMENTE
%       INVIAVEIS (o valor nao existe no catalogo) e geram viol_estrutural.
%       Consequencia: para N que nao seja potencia de 2, a fracao
%       (2^B - N)/2^B dos codigos por variavel e desperdicada.
%
%   modo_discreto = 'proporcional'  (heuristica legada, NAO consta de [DF2011])
%       indice = floor(x * N / 2^B) + 1. Elimina 100% do desperdicio, ao
%       custo de um vies estatistico (opcoes com probabilidades desiguais).
%       Mantido para comparacao experimental e testes de vies.
%
% [D2] VARIAVEL INTEIRA
%   O valor decodificado e usado DIRETAMENTE (Eq. 6), sem saturacao interna.
%   Em [DF2011] os limites do inteiro sao tratados como RESTRICAO do problema
%   (ver Eq. 12, trem de engrenagens: 12 <= z <= 60 com 6 bits => 0..63), e
%   nao por clamp no decodificador.
%   A versao anterior deste solver fazia min(x, max) e ignorava completamente
%   o limite inferior .min — bug corrigido aqui.
%
% -------------------------------------------------------------------------
% TRATAMENTO DE VALORES FORA DO DOMINIO
%   Sao sinalizados via viol_estrutural > 0. O chamador deve trata-los pela
%   regra de [DEB2000] (comparacao apenas por violacao), SEM avaliar a funcao
%   objetivo — coerente com [DEB2000] pag. 316: "It does not make sense to
%   compute the objective function value of an infeasible solution".
% -------------------------------------------------------------------------
%
% ENTRADAS
%   vetor_hibrido : vetor 1 x total_dim com a posicao da particula
%   mapa          : saida de rid_mapear_dimensoes
%   config        : struct array de configuracao das variaveis
%   modo_discreto : (opcional) 'datta' (padrao) ou 'proporcional'
%
% SAIDAS
%   vars            : vetor 1 x n_vars com os valores de projeto decodificados
%   viol_estrutural : soma das violacoes de dominio (0 = solucao avaliavel)

if nargin < 4 || isempty(modo_discreto)
    modo_discreto = 'datta';
end

n_vars          = length(config);
vars            = zeros(1, n_vars);
viol_estrutural = 0;

for k = 1:n_vars
    info = mapa(k);

    if strcmp(info.tipo, 'R')
        % Variavel real: valor lido diretamente da dimensao real
        vars(k) = vetor_hibrido(info.inicio);

    else
        bits = round(vetor_hibrido(info.inicio : info.fim));
        n    = length(bits);

        % --- Eq. (6) [DF2011]: pesos 2^(B-i), i = 1..B ---
        pesos = 2 .^ ((n-1):-1:0);
        x_int = bits(:)' * pesos(:);

        if strcmp(info.tipo, 'I')
            % ---------- [D2] INTEIRO ----------
            vars(k) = x_int;
            if x_int < config(k).min
                viol_estrutural = viol_estrutural + (config(k).min - x_int);
            elseif x_int > config(k).max
                viol_estrutural = viol_estrutural + (x_int - config(k).max);
            end

        else
            % ---------- [D1] DISCRETO ----------
            N = length(config(k).opcoes);

            if strcmp(modo_discreto, 'datta')
                % Fiel a [DF2011] Sec. 4.1: indice = x + 1, limites 1..N
                idx = x_int + 1;
                if idx > N
                    viol_estrutural = viol_estrutural + (idx - N);
                    vars(k) = config(k).opcoes(N);   % valor de preenchimento
                else
                    vars(k) = config(k).opcoes(idx);
                end

            elseif strcmp(modo_discreto, 'proporcional')
                % Heuristica legada (nao consta de [DF2011])
                n_estados = 2 ^ n;
                idx = floor(x_int * N / n_estados) + 1;
                idx = min(max(idx, 1), N);
                vars(k) = config(k).opcoes(idx);

            else
                error('rid_decodificar:modoInvalido', ...
                    'modo_discreto deve ser ''datta'' ou ''proporcional''. Recebido: ''%s''.', ...
                    modo_discreto);
            end
        end
    end
end

end
