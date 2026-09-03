function resultado = main_awruch_discreto(seed, n_runs, n_workers)
arguments
    % seed aceita tambem char/string por causa da chamada-sentinela
    % main_awruch_discreto('caso') — ver mais abaixo.
    seed {mustBeA(seed, ["double","char","string"])} = 42
    n_runs    (1,1) double {mustBePositive, mustBeInteger} = 5
    n_workers (1,1) double {mustBeNonnegative, mustBeInteger} = 0
end
% MAIN_AWRUCH_DISCRETO  Otimizacao com catalogos discretos independentes por barra.
%
% BLOCO 3 do projeto CPIO III — ORQUESTRADOR.
%
% -------------------------------------------------------------------------
% PROPOSITO E LIMITACAO
% -------------------------------------------------------------------------
% Este experimento NAO e um benchmark de validacao: o catalogo de Awruch nao
% possui artigo-fonte (ver a nota de origem na funcao local
% caso_awruch_10barras, no fim deste arquivo) e portanto NAO ha solucao
% publicada para comparar.
%
% O que ele testa e a FLEXIBILIDADE do PSO-RID: cada barra escolhe de um
% catalogo proprio, com CARDINALIDADES DIFERENTES entre si (4, 5 e 6 opcoes).
% Isso exercita o dimensionamento por variavel da Eq. (5) de [DF2011]
% (B_k^(D) = ceil(log2(N_k)), diferente para cada k), caminho que o benchmark
% de Hadi — com catalogo unico compartilhado — nao percorre.
%
% Blocos acoplados:
%   Bloco 1 (otimizador) : pso_rid.m                <- 01_pso_rid/
%   Bloco 2  (solver)    : fem_nao_linear_solver.m  <- 02_fem_nao_linear/
%
% A definicao do problema esta EMBUTIDA aqui, na funcao local
% caso_awruch_10barras. Ela herda geometria, material, apoios e cargas do
% caso Hadi via main_hadi_nao_linear('caso') — fonte unica dessa geometria
% no projeto — e troca apenas os catalogos de secoes.
%
% ACESSOR: caso = main_awruch_discreto('caso') devolve o struct do problema
% sem rodar otimizacao (usado por 06_testes/TestIntegracao.m).
%
% USO
%   main_awruch_discreto
%   main_awruch_discreto(123, 10)
%   main_awruch_discreto(123, 10, 3)   % ate 3 workers em paralelo
%
% n_workers (padrao 0 = serial) paraleliza o laco multi-start via parfor.
% Efeito colateral aceito: cada run passa a ter semente PROPRIA em vez de
% todas consumirem o mesmo stream em sequencia — ver a nota completa em
% main_hadi_nao_linear.m. preparar_pool.m aplica um teto por memoria
% disponivel, nao por nucleos.
%
% See also pso_rid, fem_nao_linear_solver, main_hadi_nao_linear

% -------------------------------------------------------------------------
% ACESSOR DO CASO (ver cabecalho)
% -------------------------------------------------------------------------
if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    garantir_caminhos({'pso_rid', 'fem_nao_linear_solver'});
    resultado = caso_awruch_10barras();
    return;
end

garantir_caminhos({'pso_rid', 'fem_nao_linear_solver'});
n_workers = preparar_pool(n_workers, n_runs);

caso = caso_awruch_10barras();

funcao_objetivo = @(areas) avaliar_projeto(areas, caso);

pso_params = struct();
pso_params.n_particulas           = 100;
pso_params.max_iter               = 1000;
pso_params.w                      = 1.0;
pso_params.c1                     = 1.0;
pso_params.c2                     = 2.0;
pso_params.pm                     = 0.15;
pso_params.auto_adaptativo        = true;
pso_params.decodificacao_discreta = 'proporcional';  % ver [D1] em pso_rid.m
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

historicos      = cell(n_runs, 1);
pesos_por_run   = nan(n_runs, 1);
areas_por_run   = cell(n_runs, 1);
viaveis_por_run = false(n_runs, 1);
linhas_log      = cell(n_runs, 1);

% Cada run e independente; parfor-safe (n_workers=0 -> serial, o padrao)
% exige semente PROPRIA por run — mesma razao/solucao de pso_rid.m (Secao C)
% e main_hadi_nao_linear.m.
parfor (r = 1:n_runs, n_workers)
    rng(seed + r - 1, 'twister');
    [areas_r, peso_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}      = hist_r;
    pesos_por_run(r)   = peso_r;
    areas_por_run{r}   = areas_r;
    viaveis_por_run(r) = (det_r.gbest_viol <= 0);

    linhas_log{r} = sprintf(['--- Run %d/%d ---\n    -> peso = %.2f kg | ' ...
        'violacao = %.3e | iters = %d\n\n'], ...
        r, n_runs, peso_r, det_r.gbest_viol, det_r.iter_executadas);
end

fprintf('%s', linhas_log{:});

melhor_peso  = inf;
melhor_areas = [];
for r = 1:n_runs
    if viaveis_por_run(r) && pesos_por_run(r) < melhor_peso
        melhor_peso  = pesos_por_run(r);
        melhor_areas = areas_por_run{r};
    end
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


% -------------------------------------------------------------------------
% >>> LOCAL caso_awruch_10barras — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_awruch_discreto('caso')
% -------------------------------------------------------------------------
function caso = caso_awruch_10barras()
% CASO_AWRUCH_10BARRAS  Trelica de 10 barras com catalogos independentes por barra.
%
% ARQUIVO DE PARAMETROS (sem logica de solver).
%
% -------------------------------------------------------------------------
% ORIGEM E STATUS DE VALIDACAO
% -------------------------------------------------------------------------
% ATENCAO: este catalogo NAO possui artigo-fonte. Foi montado pelo professor
% da disciplina apenas para EXERCITAR o solver com catalogos distintos por
% barra (em vez de um catalogo unico compartilhado, como no benchmark de
% Hadi & Alvani). Portanto:
%
%   - NAO ha solucao de referencia publicada para comparar.
%   - NAO deve ser usado como evidencia de validacao do solver.
%   - Serve como teste de FLEXIBILIDADE: confirma que o PSO-RID lida com
%     variaveis discretas de catalogos e cardinalidades diferentes entre si.
%
% CORRECAO APLICADA (Bloco 2/3): a versao anterior deste caso
% (05_legado/pre_bloco3/catalogo_awruch.m) declarava os campos
%   ref_areas = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65]
%   ref_peso  = 2325.2
% COPIADOS do benchmark de Hadi. Esses valores estavam ERRADOS neste
% contexto: nenhuma daquelas areas pertence aos catalogos definidos abaixo
% (verificado numericamente). Eram codigo morto e enganoso, e foram
% removidos. O peso de referencia fica como NaN, sinalizando "sem referencia".
%
% -------------------------------------------------------------------------
% GEOMETRIA E CARREGAMENTO
% -------------------------------------------------------------------------
% Identicos aos do benchmark de Hadi & Alvani (2003) — o que muda e apenas o
% conjunto de secoes disponiveis para cada barra. A geometria e reutilizada
% do caso Hadi, via acessor, para nao duplicar dados.
%
% SAIDA
%   caso : struct pronto para os solvers do Bloco 2, com .config_vars
%          definindo 10 variaveis discretas de catalogos independentes.

% Reaproveita geometria, material, apoios, cargas e limites do benchmark.
% main_hadi_nao_linear('caso') e a fonte unica dessa geometria no projeto.
caso = main_hadi_nao_linear('caso');

caso.nome = 'Catalogo Awruch - Trelica 10 Barras (catalogos por barra)';

% -------------------------------------------------------------------------
% CATALOGOS POR GRUPO DE BARRAS
% -------------------------------------------------------------------------
% Valores originais fornecidos em polegadas quadradas, convertidos para mm^2.
POL2_PARA_MM2 = 645.16;   % 1 in^2 = 645.16 mm^2 (exato, por definicao da polegada)

cat_A1 = [21.5, 22.5, 23.5, 24.5]             * POL2_PARA_MM2;
cat_A2 = [0.1,  0.15, 0.2,  0.25]             * POL2_PARA_MM2;
cat_A3 = [22.4, 25.4, 27.4, 29.4]             * POL2_PARA_MM2;
cat_A4 = [14.1, 14.2, 14.3, 14.4, 14.5]       * POL2_PARA_MM2;
cat_A6 = [0.5,  1.0,  1.5,  2.0,  2.5]        * POL2_PARA_MM2;
cat_A7 = [11.0, 11.3, 11.7, 12.0, 12.3, 12.5] * POL2_PARA_MM2;

caso.catalogos_por_barra = { ...
    cat_A1, cat_A2, cat_A3, cat_A4, cat_A2, ...
    cat_A6, cat_A7, cat_A3, cat_A4, cat_A2 };

% Remove os campos herdados que NAO se aplicam a este caso
caso = rmfield(caso, 'catalogo');    % nao ha catalogo unico compartilhado

% -------------------------------------------------------------------------
% SEM SOLUCAO DE REFERENCIA (ver nota de origem acima)
% -------------------------------------------------------------------------
caso.ref_areas = [];
caso.ref_peso  = NaN;

% -------------------------------------------------------------------------
% VARIAVEIS DE PROJETO (para o PSO-RID do Bloco 1)
% -------------------------------------------------------------------------
caso.config_vars = struct('tipo', {}, 'opcoes', {});
for i = 1:numel(caso.catalogos_por_barra)
    caso.config_vars(i).tipo   = 'D';
    caso.config_vars(i).opcoes = caso.catalogos_por_barra{i};
end

end


function [custo, violacao] = avaliar_projeto(areas, caso)
% [HA2003] Eq. (1) e (2); violacao separada do custo para uso com [DEB2000].
[peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);

violacao = 0;

% [HA2003] Eq. (9) — restricoes normalizadas; ver a justificativa completa
% em avaliar_projeto de main_hadi_nao_linear.m.
g_sigma = abs(Sigma(:))   / caso.sigma_max - 1;
g_desl  = abs(u_livre(:)) / caso.d_max     - 1;

violacao = violacao + sum(max(0, g_sigma));
violacao = violacao + sum(max(0, g_desl));

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


% garantir_caminhos e preparar_pool sao helpers compartilhados em
% 03_orquestrador/auxiliares/ — nao ha mais funcao local aqui.
