function resultado = main_hadi_linear(seed, n_runs)
% MAIN_HADI_LINEAR  Otimizacao da trelica de 10 barras com analise LINEAR.
%
% BLOCO 3 do projeto CPIO III — ORQUESTRADOR.
%
% Contraparte de main_hadi_nao_linear. Serve como termo de comparacao: a
% diferenca entre os dois resultados quantifica o EFEITO DA NAO LINEARIDADE
% GEOMETRICA no projeto otimo, que e a motivacao central de [HA2003].
%
% [HA2003] Sec. 6: "In order to bench mark the solutions obtained using
% non-linear analysis, an optimisation technique based on linear analysis is
% also developed. Results from the optimisation technique based on
% non-linear analysis are compared with the linear analysis."
%
% Blocos acoplados:
%   Bloco 1 (otimizador) : pso_rid.m
%   Bloco 2a (solver)    : fem_linear_solver.m
%   Bloco 2b (problema)  : problema_hadi_10barras.m
%
% USO
%   main_hadi_linear              % semente 42, 5 execucoes
%   main_hadi_linear(123, 10)
%   r = main_hadi_linear;

if nargin < 1 || isempty(seed),   seed   = 42; end
if nargin < 2 || isempty(n_runs), n_runs = 5;  end

garantir_caminhos();
rng(seed);

% 1. PROBLEMA (Bloco 2b)
caso = problema_hadi_10barras();

% 2. FUNCAO OBJETIVO (local, no fim do arquivo)
funcao_objetivo = @(areas) avaliar_projeto(areas, caso, 'linear');

% 3. CONFIGURACAO DO PSO-RID ([DF2011] Sec. 5)
pso_params = struct();
pso_params.n_particulas           = 100;
pso_params.max_iter               = 1000;
pso_params.w                      = 1.0;
pso_params.c1                     = 1.0;
pso_params.c2                     = 2.0;
pso_params.pm                     = 0.15;
pso_params.auto_adaptativo        = true;
pso_params.decodificacao_discreta = 'datta';
pso_params.tol_estagnacao         = 200;
pso_params.verbose                = true;
pso_params.print_interval         = 100;

% 4. EXECUCAO MULTI-START
imprimir_cabecalho(caso, seed, n_runs, pso_params, 'LINEAR');

historicos    = cell(n_runs, 1);
pesos_por_run = nan(n_runs, 1);
melhor_peso   = inf;
melhor_areas  = [];

for r = 1:n_runs
    fprintf('--- Executando Run %d/%d ---\n', r, n_runs);
    [areas_r, peso_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}    = hist_r;
    pesos_por_run(r) = peso_r;

    if det_r.gbest_viol <= 0 && peso_r < melhor_peso
        melhor_peso  = peso_r;
        melhor_areas = areas_r;
    end

    fprintf('    -> Run %d: peso = %.2f kg | violacao = %.3e | iters = %d | avaliacoes FEM = %d\n\n', ...
            r, peso_r, det_r.gbest_viol, det_r.iter_executadas, det_r.n_avaliacoes);
end

if isempty(melhor_areas)
    error('main_hadi_linear:semSolucaoViavel', ...
        'Nenhuma das %d execucoes encontrou solucao viavel.', n_runs);
end

% 5. AVALIACAO DETALHADA
[peso_final, Sigma_final, u_final] = fem_linear_solver(caso, melhor_areas);
[~, viol_final] = avaliar_projeto(melhor_areas, caso, 'linear');

% 6. RELATORIO E GRAFICOS
relatorio_comparativo('Trelica 10 Barras Hadi (Linear)', ...
                      melhor_areas, peso_final, ...
                      caso.ref_areas, caso.ref_peso, ...
                      Sigma_final, caso.sigma_max, ...
                      u_final, caso.d_max, viol_final);

imprimir_estatisticas_runs(pesos_por_run, caso.ref_peso);

% Confronto linear x nao linear PARA A MESMA SOLUCAO:
% mostra o quanto a analise linear subestima/superestima a resposta do
% projeto que ela mesma escolheu.
[~, Sigma_nl, u_nl] = fem_nao_linear_solver(caso, melhor_areas);
imprimir_confronto_linear_nao_linear(Sigma_final, u_final, Sigma_nl, u_nl, caso);

fig = plot_convergencia(historicos, caso.ref_peso, ...
                        'PSO-RID — Trelica 10 Barras (Analise Linear)');
salvar_figura(fig, 'convergencia_hadi_linear');

% 7. SAIDA
resultado.caso          = caso;
resultado.melhor_areas  = melhor_areas;
resultado.melhor_peso   = peso_final;
resultado.violacao      = viol_final;
resultado.Sigma         = Sigma_final;
resultado.u             = u_final;
resultado.historicos    = historicos;
resultado.pesos_por_run = pesos_por_run;
resultado.seed          = seed;

end


% #########################################################################
% FUNCOES LOCAIS
% #########################################################################

function [custo, violacao] = avaliar_projeto(areas, caso, tipo_analise)
% Ver documentacao completa em main_hadi_nao_linear.m (mesma formulacao):
% [HA2003] Eq. (1) objetivo e Eq. (2) restricoes; violacao devolvida
% separadamente para uso com a regra de [DEB2000].

if strcmp(tipo_analise, 'linear')
    [peso, Sigma, u_livre] = fem_linear_solver(caso, areas);
else
    [peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);
end

violacao = 0;

excesso_sigma = abs(Sigma) - caso.sigma_max;
excesso_sigma = excesso_sigma(excesso_sigma > 0);
violacao = violacao + sum(excesso_sigma .^ 2);

excesso_desl = abs(u_livre) - caso.d_max;
excesso_desl = excesso_desl(excesso_desl > 0);
violacao = violacao + sum(excesso_desl .^ 2);

custo = peso;
end


function imprimir_cabecalho(caso, seed, n_runs, pso_params, tipo)
fprintf('\n============================================================\n');
fprintf(' OTIMIZACAO PSO-RID — %s\n', caso.nome);
fprintf(' Analise estrutural: %s\n', tipo);
fprintf('------------------------------------------------------------\n');
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d\n', n_runs);
fprintf(' Particulas         : %d\n', pso_params.n_particulas);
fprintf(' Referencia (artigo): %.2f kg\n', caso.ref_peso);
fprintf('============================================================\n\n');
end


function imprimir_estatisticas_runs(pesos, ref_peso)
validos = pesos(isfinite(pesos));
fprintf(' ESTATISTICAS DAS EXECUCOES\n');
fprintf(' ------------------------------------------------------------\n');
fprintf('  Melhor        : %.2f kg\n', min(validos));
fprintf('  Media         : %.2f kg\n', mean(validos));
fprintf('  Desvio padrao : %.2f kg\n', std(validos));
if ~isempty(ref_peso) && ~isnan(ref_peso)
    fprintf('  Gap do melhor : %+.2f%%\n', 100*(min(validos)-ref_peso)/ref_peso);
end
fprintf(' ------------------------------------------------------------\n\n');
end


function imprimir_confronto_linear_nao_linear(S_li, u_li, S_nl, u_nl, caso)
% Reavalia o projeto otimo obtido pela analise LINEAR usando a analise NAO
% LINEAR. Se a solucao escolhida pelo modelo linear violar restricoes quando
% analisada de forma nao linear, isso demonstra na pratica o argumento
% central de [HA2003] para incluir a nao linearidade no projeto.

fprintf(' CONFRONTO: PROJETO OTIMO LINEAR REAVALIADO EM ANALISE NAO LINEAR\n');
fprintf(' ------------------------------------------------------------\n');
fprintf('  max |sigma| linear     : %9.4f MPa (limite %.2f)\n', max(abs(S_li)), caso.sigma_max);
fprintf('  max |sigma| nao linear : %9.4f MPa (limite %.2f)\n', max(abs(S_nl)), caso.sigma_max);
fprintf('  max |u| linear         : %9.4f mm  (limite %.2f)\n', max(abs(u_li)), caso.d_max);
fprintf('  max |u| nao linear     : %9.4f mm  (limite %.2f)\n', max(abs(u_nl)), caso.d_max);

viola_nl = (max(abs(S_nl)) > caso.sigma_max) || (max(abs(u_nl)) > caso.d_max);
if viola_nl
    fprintf('\n  >> O projeto otimo LINEAR VIOLA as restricoes quando reavaliado\n');
    fprintf('     em analise nao linear. Evidencia a favor de projetar\n');
    fprintf('     diretamente com o modelo nao linear ([HA2003], Sec. 1).\n');
else
    fprintf('\n  >> O projeto otimo linear permanece viavel na analise nao linear.\n');
end
fprintf(' ------------------------------------------------------------\n\n');
end


function garantir_caminhos()
if exist('pso_rid', 'file') == 2 && exist('fem_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
