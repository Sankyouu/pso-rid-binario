function resultado = main_hadi_nao_linear(seed, n_runs)
% MAIN_HADI_NAO_LINEAR  Otimizacao da trelica de 10 barras com analise NAO LINEAR.
%
% BLOCO 3 do projeto CPIO III — ORQUESTRADOR.
%
% =========================================================================
% O QUE ESTE ARQUIVO FAZ
% =========================================================================
% Este e o unico arquivo que precisa ser lido para entender o experimento
% completo. Ele amarra os tres blocos:
%
%   Bloco 1 (otimizador) : pso_rid.m                    <- 01_pso_rid/
%   Bloco 2a (solver)    : fem_nao_linear_solver.m      <- 02_fem_nao_linear/solver/
%   Bloco 2b (problema)  : problema_hadi_10barras.m     <- 02_fem_nao_linear/problemas/
%
% Toda a logica especifica do experimento (funcao objetivo, laco multi-start,
% selecao da melhor solucao) esta EMBUTIDA aqui como funcoes locais no fim do
% arquivo. Ficam fora apenas os blocos longos de impressao e grafico, que
% poluiriam a leitura do fluxo principal:
%
%   03_orquestrador/auxiliares/relatorio_comparativo.m
%   03_orquestrador/auxiliares/plot_convergencia.m
%   03_orquestrador/auxiliares/salvar_figura.m
%
% =========================================================================
% REFERENCIAS
% =========================================================================
%   [HA2003]  Hadi & Alvani (2003), Civil-Comp Press, Paper 37.
%             Benchmark de 10 barras, Sec. 6.1 e Tabela 1 (Case 2, discreto).
%   [DF2011]  Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%             Formulacao do PSO real-inteiro-discreto.
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%             Tratamento de restricoes sem parametro de penalizacao.
%
% =========================================================================
% USO
% =========================================================================
%   main_hadi_nao_linear              % semente 42, 5 execucoes
%   main_hadi_nao_linear(123)         % outra semente
%   main_hadi_nao_linear(42, 10)      % 10 execucoes
%   r = main_hadi_nao_linear;         % devolve struct com os resultados
%
% SAIDA
%   resultado : struct com .melhor_areas .melhor_peso .violacao .Sigma .u
%               .historicos .pesos_por_run .caso

if nargin < 1 || isempty(seed),   seed   = 42; end
if nargin < 2 || isempty(n_runs), n_runs = 5;  end

garantir_caminhos();
rng(seed);

% =========================================================================
% 1. PROBLEMA (Bloco 2b) — apenas parametros, nenhuma logica
% =========================================================================
caso = problema_hadi_10barras();

% =========================================================================
% 2. FUNCAO OBJETIVO
% -------------------------------------------------------------------------
% Acopla o solver FEM (Bloco 2a) as restricoes do problema. Definida como
% funcao local no fim deste arquivo (avaliar_projeto), conforme a diretriz
% de manter o orquestrador autocontido.
% =========================================================================
funcao_objetivo = @(areas) avaliar_projeto(areas, caso, 'nao_linear');

% =========================================================================
% 3. CONFIGURACAO DO PSO-RID (Bloco 1)
% -------------------------------------------------------------------------
% Valores iniciais conforme [DF2011] Sec. 5, pag. 3628. Com
% auto_adaptativo = true, cada parametro e sorteado a cada iteracao no
% intervalo [0, valor_inicial].
% =========================================================================
pso_params = struct();
pso_params.n_particulas           = 100;
pso_params.max_iter               = 1000;
pso_params.w                      = 1.0;      % [DF2011] Sec. 5
pso_params.c1                     = 1.0;      % [DF2011] Sec. 5
pso_params.c2                     = 2.0;      % [DF2011] Sec. 5
pso_params.pm                     = 0.15;     % [DF2011] Sec. 4.3/5
pso_params.auto_adaptativo        = true;     % [DF2011] Sec. 5
pso_params.decodificacao_discreta = 'datta';  % [DF2011] Sec. 4.1 (ver [D1])
pso_params.tol_estagnacao         = 200;
pso_params.verbose                = true;
pso_params.print_interval         = 100;

% =========================================================================
% 4. EXECUCAO MULTI-START
% =========================================================================
imprimir_cabecalho(caso, seed, n_runs, pso_params, 'NAO LINEAR');

historicos    = cell(n_runs, 1);
pesos_por_run = nan(n_runs, 1);
melhor_peso   = inf;
melhor_areas  = [];

for r = 1:n_runs
    fprintf('--- Executando Run %d/%d ---\n', r, n_runs);

    [areas_r, peso_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}    = hist_r;
    pesos_por_run(r) = peso_r;

    % Seleciona o melhor apenas entre solucoes VIAVEIS ([DEB2000] criterio 1)
    if det_r.gbest_viol <= 0 && peso_r < melhor_peso
        melhor_peso  = peso_r;
        melhor_areas = areas_r;
    end

    fprintf('    -> Run %d: peso = %.2f kg | violacao = %.3e | iters = %d | avaliacoes FEM = %d\n\n', ...
            r, peso_r, det_r.gbest_viol, det_r.iter_executadas, det_r.n_avaliacoes);
end

if isempty(melhor_areas)
    error('main_hadi_nao_linear:semSolucaoViavel', ...
        ['Nenhuma das %d execucoes encontrou solucao VIAVEL. ' ...
         'Aumente max_iter/n_particulas ou revise as restricoes.'], n_runs);
end

% =========================================================================
% 5. AVALIACAO DETALHADA DA MELHOR SOLUCAO
% =========================================================================
[peso_final, Sigma_final, u_final] = fem_nao_linear_solver(caso, melhor_areas);
[~, viol_final] = avaliar_projeto(melhor_areas, caso, 'nao_linear');

% =========================================================================
% 6. RELATORIO E GRAFICOS (auxiliares — blocos longos de saida)
% =========================================================================
relatorio_comparativo('Trelica 10 Barras Hadi (Nao Linear)', ...
                      melhor_areas, peso_final, ...
                      caso.ref_areas, caso.ref_peso, ...
                      Sigma_final, caso.sigma_max, ...
                      u_final, caso.d_max, viol_final);

imprimir_estatisticas_runs(pesos_por_run, caso.ref_peso);

fig = plot_convergencia(historicos, caso.ref_peso, ...
                        'PSO-RID — Trelica 10 Barras (Analise Nao Linear)');
salvar_figura(fig, 'convergencia_hadi_nao_linear');

% =========================================================================
% 7. SAIDA ESTRUTURADA
% =========================================================================
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
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################

function [custo, violacao] = avaliar_projeto(areas, caso, tipo_analise)
% AVALIAR_PROJETO  Funcao objetivo: peso sujeito a tensao e deslocamento.
%
% Formulacao do problema — [HA2003] Eq. (1) e (2):
%
%   Minimizar   W(A) = sum_i  A_i * L_i * rho                        Eq. (1)
%   sujeito a   |sigma_i| <= sigma_admissivel      i = 1..n          Eq. (2)
%               |d_j|     <= d_admissivel          j = 1..m
%
% A violacao e devolvida SEPARADAMENTE do custo (nao ha penalizacao somada
% ao peso), porque o PSO-RID trata restricoes pela regra de [DEB2000], que
% dispensa parametro de penalizacao. Ver rid_domina_deb.m.
%
% Medida de violacao: soma dos quadrados dos excessos. O quadrado torna a
% medida suave e penaliza mais fortemente violacoes grandes, ajudando a
% ordenacao entre solucoes inviaveis ([DEB2000] criterio 3).

if strcmp(tipo_analise, 'linear')
    [peso, Sigma, u_livre] = fem_linear_solver(caso, areas);
else
    [peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);
end

violacao = 0;

% Restricao de tensao (tracao e compressao)
excesso_sigma = abs(Sigma) - caso.sigma_max;
excesso_sigma = excesso_sigma(excesso_sigma > 0);
violacao = violacao + sum(excesso_sigma .^ 2);

% Restricao de deslocamento nodal
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
fprintf(' Iteracoes maximas  : %d\n', pso_params.max_iter);
fprintf(' Decodificacao      : %s\n', pso_params.decodificacao_discreta);
fprintf(' Auto-adaptativo    : %d\n', pso_params.auto_adaptativo);
fprintf(' Referencia (artigo): %.2f kg\n', caso.ref_peso);
fprintf('============================================================\n\n');
end


function imprimir_estatisticas_runs(pesos, ref_peso)
validos = pesos(isfinite(pesos));
fprintf(' ESTATISTICAS DAS EXECUCOES\n');
fprintf(' ------------------------------------------------------------\n');
fprintf('  Execucoes         : %d\n', numel(pesos));
fprintf('  Melhor            : %.2f kg\n', min(validos));
fprintf('  Media             : %.2f kg\n', mean(validos));
fprintf('  Desvio padrao     : %.2f kg\n', std(validos));
fprintf('  Pior              : %.2f kg\n', max(validos));
if ~isempty(ref_peso) && ~isnan(ref_peso)
    fprintf('  Gap do melhor     : %+.2f%% em relacao a referencia\n', ...
            100*(min(validos)-ref_peso)/ref_peso);
end
fprintf(' ------------------------------------------------------------\n\n');
end


function garantir_caminhos()
% Adiciona as pastas do projeto ao path, caso o usuario chame o orquestrador
% diretamente sem ter rodado setup_paths antes.
if exist('pso_rid', 'file') == 2 && exist('fem_nao_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
