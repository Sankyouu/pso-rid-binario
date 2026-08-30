function [best_sol, best_cost, cost_history] = pso_classico(funcao_custo, n_vars, lb, ub, pso_params)
% PSO_CLASSICO Solver PSO padrão contínuo
%
% Entradas:
%   funcao_custo : handle da função objetivo @(x) custo
%   n_vars       : número de variáveis contínuas
%   lb           : limite inferior (escalar ou vetor 1 x n_vars)
%   ub           : limite superior (escalar ou vetor 1 x n_vars)
%   pso_params   : (opcional) struct com hiperparâmetros

if nargin < 5 || isempty(pso_params)
    pso_params = struct();
end

n_particulas   = get_param(pso_params, 'n_particulas', 50);
max_iter       = get_param(pso_params, 'max_iter', 500);
w              = get_param(pso_params, 'w', 0.8);
c1             = get_param(pso_params, 'c1', 2.0);
c2             = get_param(pso_params, 'c2', 2.0);
verbose        = get_param(pso_params, 'verbose', false);

if isscalar(lb), lb = repmat(lb, 1, n_vars); end
if isscalar(ub), ub = repmat(ub, 1, n_vars); end

% Inicialização
pos = repmat(lb, n_particulas, 1) + rand(n_particulas, n_vars) .* repmat(ub - lb, n_particulas, 1);
vel = zeros(n_particulas, n_vars);
vmax = (ub - lb) / 2;

pbest_pos   = pos;
pbest_custo = inf(n_particulas, 1);

gbest_pos   = zeros(1, n_vars);
gbest_custo = inf;
cost_history = zeros(max_iter, 1);

for iter = 1:max_iter
    for i = 1:n_particulas
        custo = funcao_custo(pos(i,:));
        if custo < pbest_custo(i)
            pbest_custo(i) = custo;
            pbest_pos(i,:) = pos(i,:);
        end
        if custo < gbest_custo
            gbest_custo = custo;
            gbest_pos   = pos(i,:);
        end
    end
    cost_history(iter) = gbest_custo;

    if verbose && mod(iter, 50) == 0
        fprintf('PSO Clássico - Iter %d | Melhor Custo: %.6e\n', iter, gbest_custo);
    end

    % Atualização
    for i = 1:n_particulas
        r1 = rand(1, n_vars);
        r2 = rand(1, n_vars);
        vel(i,:) = w * vel(i,:) + c1 * r1 .* (pbest_pos(i,:) - pos(i,:)) + c2 * r2 .* (gbest_pos - pos(i,:));
        vel(i,:) = max(-vmax, min(vmax, vel(i,:)));
        pos(i,:) = pos(i,:) + vel(i,:);
        pos(i,:) = max(lb, min(ub, pos(i,:)));
    end
end

best_sol  = gbest_pos;
best_cost = gbest_custo;
end

function val = get_param(s, field_name, default_val)
if isfield(s, field_name) && ~isempty(s.(field_name))
    val = s.(field_name);
else
    val = default_val;
end
end
