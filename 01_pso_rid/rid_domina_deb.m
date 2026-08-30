function melhor = rid_domina_deb(custo_novo, viol_novo, custo_ref, viol_ref)
% RID_DOMINA_DEB  Regra de viabilidade (dominancia) de Deb, sem parametro de penalizacao.
%
% BLOCO 1 (PSO-RID) — auxiliar isolado e testavel.
%
% REFERENCIA:
%   [DEB2000] Deb, K. "An efficient constraint handling method for genetic
%             algorithms". Comput. Methods Appl. Mech. Engrg. 186 (2000)
%             311-338.  ->  pag. 316
%
% Citacao literal dos tres criterios ([DEB2000], pag. 316):
%   "1. Any feasible solution is preferred to any infeasible solution.
%    2. Among two feasible solutions, the one having better objective
%       function value is preferred.
%    3. Among two infeasible solutions, the one having smaller constraint
%       violation is preferred."
%
% Nao utiliza parametro de penalizacao (abordagem "penalty-parameter-less"),
% adotada tambem por [DF2011] Sec. 5: "the penalty-parameter-less constraint
% handling approach, proposed by Deb [8], is applied here for working with
% infeasible solutions."
%
% ENTRADAS
%   custo_novo, viol_novo : candidato
%   custo_ref,  viol_ref  : referencia (p-best, g-best, etc.)
%
% SAIDA
%   melhor : true se o candidato e PREFERIVEL a referencia
%
% CONVENCAO: violacao <= 0 significa solucao VIAVEL.

novo_viavel = (viol_novo <= 0);
ref_viavel  = (viol_ref  <= 0);

if novo_viavel && ~ref_viavel
    % Criterio 1: viavel domina inviavel
    melhor = true;

elseif ~novo_viavel && ref_viavel
    % Criterio 1: inviavel nunca domina viavel
    melhor = false;

elseif novo_viavel && ref_viavel
    % Criterio 2: ambos viaveis -> menor custo vence
    melhor = (custo_novo < custo_ref);

else
    % Criterio 3: ambos inviaveis -> menor violacao vence
    melhor = (viol_novo < viol_ref);
end

end
