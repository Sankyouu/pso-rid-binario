function resultado = main_awruch_discreto(seed, n_runs)
% MAIN_AWRUCH_DISCRETO  Otimizacao com catalogos discretos independentes por barra.
%
% BLOCO 3 do projeto CPIO III — ORQUESTRADOR.
%
% =========================================================================
% PROPOSITO E LIMITACAO
% =========================================================================
% Este experimento NAO e um benchmark de validacao: o catalogo de Awruch nao
% possui artigo-fonte (ver a nota de origem em problema_awruch_10barras.m) e
% portanto NAO ha solucao publicada para comparar.
%
% O que ele testa e a FLEXIBILIDADE do PSO-RID: cada barra escolhe de um
% catalogo proprio, com CARDINALIDADES DIFERENTES entre si (4, 5 e 6 opcoes).
% Isso exercita o dimensionamento por variavel da Eq. (5) de [DF2011]
% (B_k^(D) = ceil(log2(N_k)), diferente para cada k), caminho que o benchmark
% de Hadi — com catalogo unico compartilhado — nao percorre.
%
% Blocos acoplados:
%   Bloco 1 (otimizador) : pso_rid.m
%   Bloco 2a (solver)    : fem_nao_linear_solver.m
%   Bloco 2b (problema)  : problema_awruch_10barras.m
%
% USO
%   main_awruch_discreto
%   main_awruch_discreto(123, 10)

if nargin < 1 || isempty(seed),   seed   = 42; end
if nargin < 2 || isempty(n_runs), n_runs = 5;  end

garantir_caminhos();
rng(seed);

caso = problema_awruch_10barras();

funcao_objetivo = @(areas) avaliar_projeto(areas, caso);

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

fprintf('\n============================================================\n');
fprintf(' OTIMIZACAO PSO-RID — %s\n', caso.nome);
fprintf('------------------------------------------------------------\n');
fprintf(' Semente: %d | Execucoes: %d | Particulas: %d\n', ...
        seed, n_runs, pso_params.n_particulas);
fprintf(' SEM solucao de referencia publicada (catalogo sem artigo-fonte)\n');
fprintf('------------------------------------------------------------\n');
fprintf(' Cardinalidade dos catalogos por barra:\n   ');
fprintf('%d ', cellfun(@numel, caso.catalogos_por_barra));
fprintf('\n============================================================\n\n');

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

    fprintf('    -> Run %d: peso = %.2f kg | violacao = %.3e | iters = %d\n\n', ...
            r, peso_r, det_r.gbest_viol, det_r.iter_executadas);
end

if isempty(melhor_areas)
    error('main_awruch_discreto:semSolucaoViavel', ...
        'Nenhuma das %d execucoes encontrou solucao viavel.', n_runs);
end

[peso_final, Sigma_final, u_final] = fem_nao_linear_solver(caso, melhor_areas);
[~, viol_final] = avaliar_projeto(melhor_areas, caso);

% Sem ref_areas/ref_peso: o relatorio omite as colunas de comparacao.
relatorio_comparativo('Trelica 10 Barras (Catalogo Awruch)', ...
                      melhor_areas, peso_final, ...
                      [], NaN, ...
                      Sigma_final, caso.sigma_max, ...
                      u_final, caso.d_max, viol_final);

% Verificacao de coerencia: cada area escolhida deve pertencer ao catalogo
% daquela barra especifica.
verificar_areas_pertencem_aos_catalogos(melhor_areas, caso);

fig = plot_convergencia(historicos, NaN, 'PSO-RID — Trelica Catalogo Awruch');
salvar_figura(fig, 'convergencia_awruch_nao_linear');

resultado.caso          = caso;
resultado.melhor_areas  = melhor_areas;
resultado.melhor_peso   = peso_final;
resultado.violacao      = viol_final;
resultado.historicos    = historicos;
resultado.pesos_por_run = pesos_por_run;
resultado.seed          = seed;

end


% #########################################################################
% FUNCOES LOCAIS
% #########################################################################

function [custo, violacao] = avaliar_projeto(areas, caso)
% [HA2003] Eq. (1) e (2); violacao separada do custo para uso com [DEB2000].
[peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);

violacao = 0;

excesso_sigma = abs(Sigma) - caso.sigma_max;
excesso_sigma = excesso_sigma(excesso_sigma > 0);
violacao = violacao + sum(excesso_sigma .^ 2);

excesso_desl = abs(u_livre) - caso.d_max;
excesso_desl = excesso_desl(excesso_desl > 0);
violacao = violacao + sum(excesso_desl .^ 2);

custo = peso;
end


function verificar_areas_pertencem_aos_catalogos(areas, caso)
% Confere que a solucao respeita o catalogo ESPECIFICO de cada barra.
% Com catalogos de cardinalidades diferentes, um erro de indexacao no
% decodificador apareceria aqui.
fprintf(' VERIFICACAO DE COERENCIA COM OS CATALOGOS\n');
fprintf(' ------------------------------------------------------------\n');
tudo_ok = true;
for i = 1:numel(areas)
    cat_i = caso.catalogos_por_barra{i};
    pertence = any(abs(cat_i - areas(i)) < 1e-9);
    if ~pertence
        tudo_ok = false;
        fprintf('  Barra %2d: area %.2f NAO pertence ao seu catalogo (%d opcoes)\n', ...
                i, areas(i), numel(cat_i));
    end
end
if tudo_ok
    fprintf('  OK: todas as %d areas pertencem aos catalogos das respectivas barras.\n', ...
            numel(areas));
end
fprintf(' ------------------------------------------------------------\n\n');
end


function garantir_caminhos()
if exist('pso_rid', 'file') == 2 && exist('fem_nao_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
