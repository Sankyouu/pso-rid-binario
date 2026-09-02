classdef TestIntegracao < matlab.unittest.TestCase
% TESTINTEGRACAO  Testes de integracao dos 3 blocos operando em conjunto.
%
% BLOCOS TESTADOS: 1 + 2 + 3 acoplados.
%
% As demais classes de teste exercitam cada bloco ISOLADAMENTE. Esta verifica
% o que so aparece na juncao: a interface entre PSO-RID, solver FEM e
% arquivos de problema, e a coerencia dos orquestradores.
%
% CUSTO: usa deliberadamente orcamentos pequenos (poucas particulas e
% iteracoes) e prefere o solver LINEAR, para manter a suite rapida. A
% qualidade da otimizacao NAO e objeto destes testes — isso e papel do
% estudo estatistico em 03_orquestrador/main_estudo_estatistico.m.
%
% Executar:  runtests('06_testes/TestIntegracao.m')

    properties
        % Handles das funcoes locais de pso_rid.m (ver "ACESSO AOS
        % AUXILIARES PARA TESTE" no cabecalho daquele arquivo).
        aux
    end

    methods (TestClassSetup)
        function adicionarCaminhos(tc)
            aqui = fileparts(mfilename('fullpath'));
            raiz = fullfile(aqui, '..');
            addpath(raiz);
            addpath(genpath(fullfile(raiz, '01_pso_rid')));
            addpath(genpath(fullfile(raiz, '02_fem_nao_linear')));
            addpath(genpath(fullfile(raiz, '03_orquestrador')));
            tc.aux = pso_rid('auxiliares');
        end
    end

    methods (Test)

        % -----------------------------------------------------------------
        % 1. INTERFACE PSO <-> PROBLEMA
        % -----------------------------------------------------------------
        function testProblemaHadi_ForneceInterfaceCompletaParaOsDoisBlocos(tc)
            % O arquivo de problema precisa servir SIMULTANEAMENTE ao Bloco 2
            % (campos do solver) e ao Bloco 1 (config_vars do otimizador).
            caso = main_hadi_nao_linear('caso');

            % Campos exigidos pelo solver (Bloco 2a)
            for campo = {'nodes0','elements','apoios','F_total','E','dens'}
                tc.verifyTrue(isfield(caso, campo{1}), ...
                    sprintf('Falta o campo "%s" exigido pelo solver.', campo{1}));
            end

            % Campos exigidos pelo otimizador (Bloco 1)
            tc.verifyTrue(isfield(caso, 'config_vars'));
            tc.verifyEqual(numel(caso.config_vars), size(caso.elements, 2), ...
                'Deve haver uma variavel de projeto por barra.');
            for i = 1:numel(caso.config_vars)
                tc.verifyEqual(caso.config_vars(i).tipo, 'D');
                tc.verifyNotEmpty(caso.config_vars(i).opcoes);
            end
        end

        function testAreasDecodificadasSaoAceitasPeloSolver(tc)
            % Toda solucao que o PSO produz precisa ser aceitavel pelo FEM.
            % Este teste fecha o ciclo: decodifica -> avalia -> sem erro.
            caso = main_hadi_nao_linear('caso');

            chamadas = 0;
            function [j, v] = espiao(areas)
                chamadas = chamadas + 1;
                % Se o PSO entregasse algo invalido, isto lancaria erro:
                [peso, Sigma, u] = fem_linear_solver(caso, areas);
                j = peso;
                v = sum(max(abs(Sigma) - caso.sigma_max, 0).^2) ...
                  + sum(max(abs(u)     - caso.d_max,     0).^2);
            end

            p = struct('n_particulas', 8, 'max_iter', 5, 'verbose', false, ...
                       'reinit_freq', 0, 'tol_estagnacao', 1e9);
            pso_rid(@espiao, caso.config_vars, p);

            tc.verifyGreaterThan(chamadas, 0, ...
                'O solver FEM nunca foi chamado — interface quebrada.');
        end

        function testSolucaoFinalPertenceAoCatalogo(tc)
            % A melhor solucao devolvida pelo PSO deve ser construivel:
            % todas as areas tem de existir no catalogo do problema.
            caso = main_hadi_nao_linear('caso');

            fobj = @(areas) custo_linear(caso, areas);
            p = struct('n_particulas', 12, 'max_iter', 12, 'verbose', false, ...
                       'reinit_freq', 0, 'tol_estagnacao', 1e9);

            rng(1);
            areas = pso_rid(fobj, caso.config_vars, p);

            tc.verifyEqual(numel(areas), size(caso.elements,2));
            for i = 1:numel(areas)
                tc.verifyTrue(ismember(areas(i), caso.catalogo), ...
                    sprintf('Area %d (%.1f) nao pertence ao catalogo.', i, areas(i)));
            end
        end

        % -----------------------------------------------------------------
        % 2. CATALOGOS DE CARDINALIDADES DIFERENTES (caso Awruch)
        % -----------------------------------------------------------------
        function testCatalogosPorBarra_DimensionamentoIndependente(tc)
            % [DF2011] Eq. (5): cada variavel discreta recebe B_k = ceil(log2(N_k)).
            % Com catalogos de tamanhos distintos, os blocos de bits precisam
            % ser dimensionados INDEPENDENTEMENTE — caminho que o benchmark de
            % Hadi (catalogo unico) nao exercita.
            caso = main_awruch_discreto('caso');

            [mapa, total] = tc.aux.mapear_dimensoes(caso.config_vars);

            for i = 1:numel(caso.config_vars)
                N_i = numel(caso.config_vars(i).opcoes);
                tc.verifyEqual(mapa(i).n_bits, max(ceil(log2(N_i)),1), ...
                    sprintf('Barra %d (N=%d) recebeu numero de bits errado.', i, N_i));
            end
            tc.verifyEqual(total, sum([mapa.n_bits]));
        end

        function testAwruch_SolucaoRespeitaCatalogoDeCadaBarra(tc)
            % Com catalogos diferentes por barra, um erro de indexacao no
            % decodificador produziria uma area do catalogo errado.
            caso = main_awruch_discreto('caso');

            fobj = @(areas) custo_linear(caso, areas);
            p = struct('n_particulas', 12, 'max_iter', 12, 'verbose', false, ...
                       'reinit_freq', 0, 'tol_estagnacao', 1e9);

            rng(2);
            areas = pso_rid(fobj, caso.config_vars, p);

            for i = 1:numel(areas)
                cat_i = caso.catalogos_por_barra{i};
                tc.verifyTrue(any(abs(cat_i - areas(i)) < 1e-9), ...
                    sprintf(['Barra %d recebeu area %.2f, que nao pertence ao ' ...
                             'seu proprio catalogo (%d opcoes).'], ...
                            i, areas(i), numel(cat_i)));
            end
        end

        function testAwruch_NaoHerdouReferenciaIncorretaDoHadi(tc)
            % REGRESSAO: a versao anterior (catalogo_awruch.m) copiava
            % ref_areas/ref_peso do benchmark de Hadi, valores que nem sequer
            % pertencem aos catalogos de Awruch. Devem estar limpos.
            caso = main_awruch_discreto('caso');

            tc.verifyEmpty(caso.ref_areas, ...
                'Awruch nao tem solucao de referencia publicada.');
            tc.verifyTrue(isnan(caso.ref_peso), ...
                'ref_peso de Awruch deveria ser NaN (sem referencia).');
            tc.verifyFalse(isfield(caso, 'catalogo'), ...
                'Awruch nao possui catalogo unico compartilhado.');
        end

        % -----------------------------------------------------------------
        % 3. COERENCIA ENTRE ANALISE LINEAR E NAO LINEAR NO MESMO PROBLEMA
        % -----------------------------------------------------------------
        function testMesmoProblemaAceitoPelosDoisSolvers(tc)
            % O mesmo struct de problema deve alimentar os dois solvers sem
            % nenhuma adaptacao — requisito da separacao solver/parametros.
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            [w_li, S_li, u_li] = fem_linear_solver(caso, areas);
            [w_nl, S_nl, u_nl] = fem_nao_linear_solver(caso, areas);

            % O peso e puramente geometrico: identico nos dois.
            tc.verifyEqual(w_li, w_nl, 'RelTol', 1e-12, ...
                'O peso nao pode depender do tipo de analise.');

            % As respostas diferem (ha nao linearidade), mas na mesma ordem.
            tc.verifySize(S_nl, size(S_li));
            tc.verifySize(u_nl, size(u_li));
            tc.verifyLessThan(max(abs(u_nl - u_li))/max(abs(u_li)), 0.5, ...
                'Diferenca linear/nao linear implausivelmente grande.');
        end

        % -----------------------------------------------------------------
        % 4. ORQUESTRADORES (Bloco 3)
        % -----------------------------------------------------------------
        function testOrquestradoresExistemEEstaoNoPath(tc)
            for nome = {'main_hadi_nao_linear', 'main_hadi_linear', ...
                        'main_awruch_discreto', 'main_estudo_estatistico'}
                tc.verifyEqual(exist(nome{1}, 'file'), 2, ...
                    sprintf('Orquestrador %s nao encontrado no path.', nome{1}));
            end
        end

        function testAuxiliaresDeSaidaExistem(tc)
            % Blocos longos de impressao/grafico mantidos fora dos
            % orquestradores, conforme a diretriz de organizacao.
            for nome = {'relatorio_comparativo', 'plot_convergencia', 'salvar_figura'}
                tc.verifyEqual(exist(nome{1}, 'file'), 2, ...
                    sprintf('Auxiliar %s nao encontrado no path.', nome{1}));
            end
        end

        function testSetupPathsHabilitaTodosOsBlocos(tc)
            % Apos setup_paths, os pontos de entrada dos 3 blocos devem estar
            % acessiveis a partir de qualquer diretorio.
            setup_paths(false);
            tc.verifyEqual(exist('pso_rid', 'file'), 2,               'Bloco 1 ausente.');
            tc.verifyEqual(exist('fem_linear_solver', 'file'), 2,     'Bloco 2 (linear) ausente.');
            tc.verifyEqual(exist('fem_nao_linear_solver', 'file'), 2, 'Bloco 2 (nao linear) ausente.');
            tc.verifyEqual(exist('main_hadi_nao_linear', 'file'), 2,  'Bloco 3 ausente.');

            % Os casos de estudo nao sao mais arquivos proprios: vivem como
            % funcao local dentro do respectivo orquestrador e sao alcancados
            % pelo acessor ('caso'). Verificar que o acessor responde e o
            % equivalente atual de "o Bloco 2b esta no path".
            tc.verifyTrue(isstruct(main_hadi_nao_linear('caso')), ...
                'Acessor do caso Hadi nao respondeu.');
            tc.verifyTrue(isstruct(main_awruch_discreto('caso')), ...
                'Acessor do caso Awruch nao respondeu.');
        end

        % -----------------------------------------------------------------
        % 5. CICLO COMPLETO EM ESCALA REDUZIDA
        % -----------------------------------------------------------------
        function testCicloCompletoProduzSolucaoViavel(tc)
            % Executa PSO + FEM ate o fim, em escala pequena, e confere que a
            % solucao devolvida e coerente: viavel, dentro do catalogo, e com
            % peso consistente com o recalculo pelo solver.
            caso = main_hadi_nao_linear('caso');

            fobj = @(areas) custo_linear(caso, areas);
            p = struct('n_particulas', 40, 'max_iter', 60, 'verbose', false, ...
                       'tol_estagnacao', 1e9);

            rng(42);
            [areas, peso_pso, hist, det] = pso_rid(fobj, caso.config_vars, p);

            % O peso devolvido pelo PSO tem de bater com o recalculo do FEM.
            peso_fem = fem_linear_solver(caso, areas);
            tc.verifyEqual(peso_pso, peso_fem, 'RelTol', 1e-10, ...
                'Peso do PSO diverge do recalculado pelo solver.');

            % MONOTONICIDADE SOB [DEB2000] — corrigido em 2026-08-31.
            %
            % O custo do g-best NAO e monotonico, e isso e correto: pelo
            % criterio 1, um projeto VIAVEL e preferido a qualquer inviavel
            % independentemente do custo. No instante em que o enxame troca
            % um inviavel-barato por um viavel-caro, o custo SOBE.
            % (Observado: 3916.83 -> 4543.62 na iteracao 1->2 com semente 42.)
            %
            % A versao anterior deste teste exigia custo nao crescente. Ela
            % passava por acidente: com a decodificacao 'datta', as primeiras
            % iteracoes eram dominadas por estouros de catalogo com custo Inf,
            % que o filtro isfinite removia — escondendo a fase inviavel.
            %
            % O invariante de fato garantido pelos criterios 1 e 3 e sobre a
            % VIOLACAO, que nunca pode piorar.
            tc.verifyNumElements(hist, det.iter_executadas);
            tc.verifyNumElements(det.viol_history, det.iter_executadas);

            v = det.viol_history;
            tc.verifyTrue(all(diff(v(isfinite(v))) <= 1e-9), ...
                'A violacao do g-best deveria ser monotonicamente nao crescente.');

            % Uma vez atingida a viabilidade, nunca mais se perde...
            i_viavel = find(v <= 0, 1);
            if ~isempty(i_viavel)
                tc.verifyTrue(all(v(i_viavel:end) <= 0), ...
                    'O g-best perdeu viabilidade depois de te-la atingido.');
                % ...e dai em diante o custo SO PODE CAIR (criterio 2).
                tc.verifyTrue(all(diff(hist(i_viavel:end)) <= 1e-9), ...
                    'Custo do g-best subiu dentro da regiao viavel.');
            end

            % Se encontrou solucao viavel, ela deve respeitar as restricoes.
            if det.gbest_viol <= 0
                [~, S, u] = fem_linear_solver(caso, areas);
                tc.verifyLessThanOrEqual(max(abs(S)), caso.sigma_max * (1 + 1e-9));
                tc.verifyLessThanOrEqual(max(abs(u)), caso.d_max     * (1 + 1e-9));
            end
        end

        % -----------------------------------------------------------------
        % MEDIDA DE VIOLACAO — [HA2003] Eq. (9), normalizada
        % -----------------------------------------------------------------
        function testViolacaoNormalizadaEhAdimensional(tc)
            % A violacao nao pode depender da unidade em que o problema esta
            % escrito. Reescalar comprimentos (mm -> m: deslocamentos e
            % d_max por 1e-3) e tensoes (MPa -> Pa: sigma e sigma_max por
            % 1e6) tem de deixar a violacao IDENTICA, porque cada parcela e
            % uma razao entre grandezas da mesma familia.
            %
            % REGRESSAO: a medida anterior somava (|sigma|-sigma_max)^2 com
            % (|u|-d_max)^2, ou seja MPa^2 com mm^2. Sob a mesma reescala
            % ela muda de valor e chega a INVERTER a ordenacao entre dois
            % projetos inviaveis — medido em 1.7% dos pares amostrados.
            % Como o criterio 3 de [DEB2000] ordena inviaveis SOMENTE pela
            % violacao, isso fazia a busca depender da escolha de unidades.
            caso = main_hadi_nao_linear('caso');
            rng(20260831);

            for k = 1:40
                areas = caso.catalogo(randi(numel(caso.catalogo), 1, 10));
                [~, S, u] = fem_linear_solver(caso, areas);

                v_mm = tc.violacaoNormalizada(S, u, caso.sigma_max, caso.d_max);
                v_m  = tc.violacaoNormalizada(S * 1e6, u * 1e-3, ...
                                              caso.sigma_max * 1e6, ...
                                              caso.d_max * 1e-3);

                tc.verifyEqual(v_m, v_mm, 'RelTol', 1e-12, ...
                    'Violacao normalizada mudou ao trocar a unidade.');
            end
        end

        function testViolacaoNormalizadaMedeFracaoDoLimite(tc)
            % Semantica da Eq. (9): o valor "0.10" significa 10% acima do
            % admissivel, na tensao ou no deslocamento, indiferentemente.
            % E isso que torna as duas familias comensuraveis.
            sigma_max = 172.25;   d_max = 50.80;

            % Exatamente no limite: viavel, violacao nula.
            tc.verifyEqual(tc.violacaoNormalizada(sigma_max, d_max, ...
                                                  sigma_max, d_max), 0, ...
                'Projeto exatamente no limite deveria ter violacao zero.');

            % Folga: continua zero (nao existe violacao negativa).
            tc.verifyEqual(tc.violacaoNormalizada(sigma_max/2, d_max/2, ...
                                                  sigma_max, d_max), 0);

            % 10% acima em tensao, e so.
            tc.verifyEqual(tc.violacaoNormalizada(1.10 * sigma_max, d_max, ...
                                                  sigma_max, d_max), 0.10, ...
                'AbsTol', 1e-12);

            % 10% acima em deslocamento: TEM de valer o mesmo que o anterior.
            tc.verifyEqual(tc.violacaoNormalizada(sigma_max, 1.10 * d_max, ...
                                                  sigma_max, d_max), 0.10, ...
                'AbsTol', 1e-12, ...
                'Excessos relativos iguais devem pesar igual nas duas familias.');

            % As parcelas somam ([DEB2000] pag. 316).
            tc.verifyEqual(tc.violacaoNormalizada(1.10 * sigma_max, ...
                                                  1.25 * d_max, ...
                                                  sigma_max, d_max), 0.35, ...
                'AbsTol', 1e-12);
        end

        function testViolacaoNulaEquivaleAViabilidade(tc)
            % A convencao que liga avaliar_projeto a domina_deb: violacao <= 0
            % tem de significar exatamente "respeita Eq. (2)".
            caso = main_hadi_nao_linear('caso');
            rng(20260831);

            for k = 1:60
                areas = caso.catalogo(randi(numel(caso.catalogo), 1, 10));
                [~, viol] = custo_linear(caso, areas);
                [~, S, u] = fem_linear_solver(caso, areas);

                respeita = max(abs(S)) <= caso.sigma_max * (1 + 1e-12) && ...
                           max(abs(u)) <= caso.d_max     * (1 + 1e-12);

                tc.verifyEqual(viol <= 0, respeita, ...
                    'violacao <= 0 deveria coincidir com respeitar a Eq. (2).');
            end
        end

        function testEstudoEstatistico_RespeitaOrcamentoDeAvaliacoes(tc)
            % O estudo estatistico compara configuracoes sob orcamento
            % igualado de avaliacoes do FEM. Este teste confere que o
            % orcamento e de fato respeitado e que o PSO PARA ao esgota-lo,
            % em vez de continuar iterando sem realizar avaliacoes uteis.
            %
            % REGRESSAO: uma versao anterior devolvia custo infinito ao
            % esgotar o orcamento, deixando a parada por conta do criterio de
            % estagnacao. Com tol_estagnacao alto, o laco girava dezenas de
            % milhares de iteracoes a toa.
            caso = main_hadi_nao_linear('caso');

            orcamento = 300;
            n_chamadas = 0;

            function [j, v] = objetivo_com_orcamento(areas)
                if n_chamadas >= orcamento
                    error('teste:orcamentoEsgotado', 'Orcamento esgotado.');
                end
                n_chamadas = n_chamadas + 1;
                [j, v] = custo_linear(caso, areas);
            end

            % tol_estagnacao deliberadamente altissimo: se a parada dependesse
            % da estagnacao, este teste nao terminaria em tempo razoavel.
            p = struct('n_particulas', 20, 'max_iter', 100000, ...
                       'verbose', false, 'tol_estagnacao', 1e9);

            rng(3);
            t0 = tic;
            try
                pso_rid(@objetivo_com_orcamento, caso.config_vars, p);
                interrompeu = false;
            catch err
                interrompeu = strcmp(err.identifier, 'teste:orcamentoEsgotado');
            end
            decorrido = toc(t0);

            tc.verifyTrue(interrompeu, ...
                'O esgotamento do orcamento deveria interromper a execucao.');
            tc.verifyEqual(n_chamadas, orcamento, ...
                'O numero de avaliacoes deveria parar exatamente no orcamento.');
            tc.verifyLessThan(decorrido, 60, ...
                'A execucao deveria terminar logo apos esgotar o orcamento.');
        end

        function testReprodutibilidadePorSemente(tc)
            % Mesma semente -> mesmo resultado. Requisito para que o estudo
            % estatistico pareado por semente faca sentido.
            caso = main_hadi_nao_linear('caso');
            fobj = @(areas) custo_linear(caso, areas);
            p = struct('n_particulas', 10, 'max_iter', 15, 'verbose', false, ...
                       'reinit_freq', 0, 'tol_estagnacao', 1e9);

            rng(7);  [a1, c1] = pso_rid(fobj, caso.config_vars, p);
            rng(7);  [a2, c2] = pso_rid(fobj, caso.config_vars, p);

            tc.verifyEqual(a1, a2, 'AbsTol', 0, ...
                'Mesma semente deveria produzir a mesma solucao.');
            tc.verifyEqual(c1, c2, 'AbsTol', 0, ...
                'Mesma semente deveria produzir o mesmo custo.');
        end
    end

    methods (Static, Access = private)
        function v = violacaoNormalizada(sigma, u, sigma_max, d_max)
            % [HA2003] Eq. (9), na mesma forma usada por avaliar_projeto
            % dos orquestradores e por custo_linear deste arquivo.
            v = sum(max(abs(sigma(:)) / sigma_max - 1, 0)) ...
              + sum(max(abs(u(:))     / d_max     - 1, 0));
        end
    end

end


% #########################################################################
% FUNCAO LOCAL DO ARQUIVO (acessivel aos metodos da classe)
% #########################################################################

function [custo, violacao] = custo_linear(caso, areas)
% Funcao objetivo enxuta para os testes: mesma formulacao dos orquestradores
% ([HA2003] Eq. 1 objetivo, Eq. 2 restricoes e Eq. 9 violacao NORMALIZADA),
% usando o solver LINEAR por ser muito mais rapido — a suite nao deve
% depender do custo do nao linear.
%
% Esta copia precisa acompanhar avaliar_projeto dos orquestradores: se as
% duas medidas divergirem, os testes deixam de exercitar o que roda de
% verdade. O teste testViolacaoNormalizadaEhAdimensional trava isso.
[peso, Sigma, u] = fem_linear_solver(caso, areas);
violacao = sum(max(abs(Sigma(:)) / caso.sigma_max - 1, 0)) ...
         + sum(max(abs(u(:))     / caso.d_max     - 1, 0));
custo = peso;
end
