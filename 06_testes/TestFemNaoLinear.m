classdef TestFemNaoLinear < matlab.unittest.TestCase
% TESTFEMNAOLINEAR  Testes do solver geometricamente nao linear (Bloco 2a).
%
% BLOCO TESTADO: 2 (FEM nao linear), isoladamente — nenhum PSO envolvido.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%   [HA2003]  Hadi & Alvani (2003), Civil-Comp Press, Paper 37.
%             Eq. (13) KT = KE + KG;  Eq. (14) matriz geometrica;
%             Eq. (18) passo de Newton-Raphson;  Sec. 6.1 e Tabela 1.
%   [CRI1991] Crisfield, Vol. 1 (1991) — formulacao de base citada por [HA2003].
%
% -------------------------------------------------------------------------
% ESTRATEGIA DE VALIDACAO (do mais forte para o mais fraco)
% -------------------------------------------------------------------------
%   1. SOLUCAO ANALITICA EXATA — trelica rasa de 2 barras. Compara o
%      deslocamento do apice com a formula fechada derivada para a MESMA
%      medida de deformacao usada pelo solver. Independente de terceiros.
%   2. EQUILIBRIO GLOBAL — residuo F_ext - F_int deve ser ~0.
%   3. LIMITE LINEAR — com carga tendendo a zero, o nao linear deve
%      convergir para o solver linear ([HA2003] Sec. 4: "initially no change
%      in geometry is registered, KG = 0, rendering this stage similar to a
%      linear analysis").
%   4. INVARIANCIA NUMERICA — o resultado nao pode depender do numero de
%      incrementos de carga.
%   5. CARACTERIZACAO DO BENCHMARK — reproduz Hadi Case 2 e DOCUMENTA as
%      discrepancias conhecidas (nao e um teste de igualdade estrita).
%
% Executar:  runtests('06_testes/TestFemNaoLinear.m')

    methods (TestClassSetup)
        function adicionarCaminhos(~)
            aqui = fileparts(mfilename('fullpath'));
            raiz = fullfile(aqui, '..');
            addpath(genpath(fullfile(raiz, '01_pso_rid')));
            addpath(fullfile(raiz, '05_legado'));
            addpath(genpath(fullfile(raiz, '02_fem_nao_linear')));
            % O caso Hadi vem do acessor main_hadi_nao_linear('caso') e o
            % oraculo analitico (problema_trelica_rasa_2barras) vive aqui
            % em 06_testes, ao lado desta classe.
            addpath(genpath(fullfile(raiz, '03_orquestrador')));
            addpath(aqui);
        end
    end

    methods (Test)

        % -----------------------------------------------------------------
        % 1. VALIDACAO CONTRA SOLUCAO ANALITICA EXATA
        % -----------------------------------------------------------------
        function testAnalitico_DeslocamentoDoApiceBateComFormulaFechada(tc)
            % Trelica rasa de 2 barras: o deslocamento do apice calculado
            % pelo Newton-Raphson deve coincidir com a solucao fechada
            %
            %     P(w) = -2 * E*A*((L(w)-L0)/L0) * (h-w)/L(w)
            %
            % derivada em problema_trelica_rasa_2barras.m para deformacao de
            % engenharia — a mesma medida usada pelo solver.
            %
            % Este e o teste mais forte da suite: valida simultaneamente a
            % atualizacao da geometria, a matriz geometrica KG (Eq. 14),
            % a montagem das forcas internas e a convergencia do N-R.
            caso = problema_trelica_rasa_2barras();

            fracoes = [0.10 0.25 0.50 0.75 0.90];
            for f = fracoes
                P       = f * caso.P_critica;
                w_exato = caso.w_de_P(P);

                prob = caso;
                prob.F_total    = zeros(6,1);
                prob.F_total(caso.dof_apice_v) = -P;

                [~, ~, ~, d] = fem_nao_linear_solver(prob, caso.areas, ...
                                    struct('n_inc', 20, 'tol', 1e-10));

                w_fem = -d.u(caso.dof_apice_v);   % positivo para baixo

                tc.verifyTrue(d.convergiu, ...
                    sprintf('N-R nao convergiu em %.0f%% da carga critica.', 100*f));
                tc.verifyEqual(w_fem, w_exato, 'RelTol', 1e-8, ...
                    sprintf(['Deslocamento do apice diverge da solucao exata ' ...
                             'em %.0f%% da carga critica.'], 100*f));
            end
        end

        function testAnalitico_RegimeFortementeNaoLinear(tc)
            % Em 90% da carga critica o deslocamento do apice ja e ~28% da
            % altura inicial: regime claramente nao linear, onde uma analise
            % linear erraria de forma grosseira. O acordo com a solucao
            % exata aqui demonstra que a nao linearidade esta sendo captada
            % (e nao apenas que ambos concordam no regime trivial).
            caso = problema_trelica_rasa_2barras();
            P    = 0.90 * caso.P_critica;

            prob = caso;
            prob.F_total = zeros(6,1);
            prob.F_total(caso.dof_apice_v) = -P;

            [~,~,~,d_nl] = fem_nao_linear_solver(prob, caso.areas, ...
                                struct('n_inc', 20, 'tol', 1e-10));
            [~,~,~,d_li] = fem_linear_solver(prob, caso.areas);

            w_nl    = -d_nl.u(caso.dof_apice_v);
            w_li    = -d_li.u(caso.dof_apice_v);
            w_exato = caso.w_de_P(P);

            % O nao linear acerta...
            tc.verifyEqual(w_nl, w_exato, 'RelTol', 1e-8);

            % ...e o linear erra substancialmente (prova que o caso e
            % genuinamente nao linear, nao um teste trivial).
            erro_linear = abs(w_li - w_exato) / w_exato;
            tc.verifyGreaterThan(erro_linear, 0.20, ...
                ['A analise linear deveria errar >20%% neste regime; se nao ' ...
                 'errar, o caso de teste nao esta exercitando a nao linearidade.']);

            fprintf(['\n[CARACTERIZACAO NAO-LINEARIDADE] P = 90%% da carga critica\n' ...
                     '  w exato        : %8.4f mm\n' ...
                     '  w nao linear   : %8.4f mm  (erro %.2e)\n' ...
                     '  w linear       : %8.4f mm  (erro %.1f%%)\n'], ...
                    w_exato, w_nl, abs(w_nl-w_exato)/w_exato, w_li, 100*erro_linear);
        end

        function testAnalitico_OraculoRejeitaCargaAcimaDaCritica(tc)
            % A trelica rasa tem ponto limite (snap-through). O oraculo deve
            % recusar explicitamente cargas sem equilibrio no ramo estavel,
            % em vez de devolver um numero silenciosamente errado.
            caso = problema_trelica_rasa_2barras();
            tc.verifyError(@() caso.w_de_P(1.01 * caso.P_critica), ...
                'problema_trelica_rasa_2barras:cargaAcimaDaCritica');
        end

        % -----------------------------------------------------------------
        % 2. EQUILIBRIO GLOBAL
        % -----------------------------------------------------------------
        function testEquilibrioGlobal_ResiduoDesprezivel(tc)
            % Ao fim da analise, as forcas internas devem equilibrar as
            % externas nos GDLs livres.
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            [~,~,~,d] = fem_nao_linear_solver(caso, areas, ...
                            struct('n_inc', 10, 'tol', 1e-10));

            F_int = tc.montarForcasInternas(caso, d);
            R = caso.F_total - F_int;
            R(caso.apoios) = 0;

            escala = norm(caso.F_total);
            tc.verifyLessThan(norm(R)/escala, 1e-8, ...
                'Residuo de equilibrio global acima da tolerancia.');
        end

        function testReacoesEquilibramCargaAplicada(tc)
            % A soma das reacoes nos apoios deve equilibrar a resultante
            % externa (equilibrio global de corpo livre).
            caso  = main_hadi_nao_linear('caso');
            [~,~,~,d] = fem_nao_linear_solver(caso, caso.ref_areas, ...
                            struct('n_inc', 10, 'tol', 1e-10));

            F_int = tc.montarForcasInternas(caso, d);

            % Nos apoios, a forca interna e a reacao.
            soma_vertical_reacoes = sum(F_int(caso.apoios(mod(caso.apoios,2)==0)));
            carga_vertical_total  = sum(caso.F_total);

            tc.verifyEqual(soma_vertical_reacoes, -carga_vertical_total, ...
                'RelTol', 1e-6, ...
                'Reacoes verticais nao equilibram a carga aplicada.');
        end

        % -----------------------------------------------------------------
        % 3. LIMITE LINEAR
        % -----------------------------------------------------------------
        function testCargaPequena_ConvergeParaSolucaoLinear(tc)
            % [HA2003] Sec. 4: com deslocamentos despreziveis, KG -> 0 e a
            % analise nao linear degenera na linear.
            %
            % FORMULACAO DO TESTE: em vez de exigir coincidencia a uma
            % tolerancia fixa arbitraria, verifica-se a propriedade que de
            % fato caracteriza a convergencia — a discrepancia entre as duas
            % analises DIMINUI proporcionalmente ao fator de carga.
            %
            % (Uma tolerancia fixa e fragil aqui: para cargas muito pequenas
            % a diferenca residual passa a ser dominada pela tolerancia de
            % convergencia do Newton-Raphson, e nao pela fisica. Testar a
            % TENDENCIA e mais significativo e mais robusto.)
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            fatores    = [1e-1, 1e-2, 1e-3];
            discrep_u  = zeros(size(fatores));

            for k = 1:numel(fatores)
                prob = caso;
                prob.F_total = caso.F_total * fatores(k);

                [~, ~, u_nl] = fem_nao_linear_solver(prob, areas, ...
                                   struct('n_inc', 5, 'tol', 1e-12));
                [~, ~, u_li] = fem_linear_solver(prob, areas);

                discrep_u(k) = max(abs(u_nl - u_li)) / max(abs(u_li));
            end

            % A discrepancia deve encolher ao reduzir a carga...
            for k = 2:numel(fatores)
                tc.verifyLessThan(discrep_u(k), discrep_u(k-1), ...
                    sprintf(['A discrepancia nao-linear/linear deveria diminuir ' ...
                             'ao reduzir a carga (fator %g -> %g).'], ...
                            fatores(k-1), fatores(k)));
            end

            % ...e ser aproximadamente proporcional ao fator de carga
            % (convergencia de primeira ordem): ao dividir a carga por 10,
            % a discrepancia deve cair por um fator proximo de 10.
            razao = discrep_u(1) / discrep_u(end);
            tc.verifyGreaterThan(razao, 50, ...
                ['Reduzir a carga em 100x deveria reduzir a discrepancia em ' ...
                 'ordem semelhante (convergencia de 1a ordem para o linear).']);

            % Em carga muito pequena, as duas analises praticamente coincidem.
            tc.verifyLessThan(discrep_u(end), 1e-4, ...
                'Com carga 1000x menor, nao linear e linear deveriam coincidir.');

            fprintf(['\n[CARACTERIZACAO LIMITE LINEAR] discrepancia relativa ' ...
                     'max|u_nl - u_li| / max|u_li|:\n' ...
                     '  fator de carga 1e-1 -> %.3e\n' ...
                     '  fator de carga 1e-2 -> %.3e\n' ...
                     '  fator de carga 1e-3 -> %.3e\n'], ...
                    discrep_u(1), discrep_u(2), discrep_u(3));
        end

        % -----------------------------------------------------------------
        % 4. INVARIANCIA NUMERICA
        % -----------------------------------------------------------------
        function testResultadoIndependeDoNumeroDeIncrementos(tc)
            % Um Newton-Raphson corretamente implementado converge para o
            % MESMO equilibrio independentemente de quantos incrementos de
            % carga sao usados (o incremento afeta a robustez, nao a solucao).
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            [~, S_ref, u_ref] = fem_nao_linear_solver(caso, areas, ...
                                    struct('n_inc', 10, 'tol', 1e-10));

            for n_inc = [1 5 20 50]
                [~, S, u] = fem_nao_linear_solver(caso, areas, ...
                                struct('n_inc', n_inc, 'tol', 1e-10));
                tc.verifyEqual(u, u_ref, 'AbsTol', 1e-6, ...
                    sprintf('Deslocamentos mudaram com n_inc = %d.', n_inc));
                tc.verifyEqual(S, S_ref, 'AbsTol', 1e-6, ...
                    sprintf('Tensoes mudaram com n_inc = %d.', n_inc));
            end
        end

        function testConvergenciaEhReportada(tc)
            % O solver deve informar se convergiu. Com tolerancia
            % inatingivel em 1 iteracao, deve reportar convergiu = false
            % em vez de devolver um resultado silenciosamente nao convergido.
            caso = main_hadi_nao_linear('caso');
            [~,~,~,d] = fem_nao_linear_solver(caso, caso.ref_areas, ...
                            struct('n_inc', 1, 'max_iter', 1, 'tol', 1e-14));
            tc.verifyFalse(d.convergiu, ...
                'Solver deveria reportar nao convergencia com max_iter = 1.');
        end

        % -----------------------------------------------------------------
        % 5. CARACTERIZACAO DO BENCHMARK DE HADI
        % -----------------------------------------------------------------
        function testHadi_CaracterizacaoDaSolucaoDeReferencia(tc)
            % TESTE DE CARACTERIZACAO — trava os valores atualmente obtidos
            % e DOCUMENTA as discrepancias conhecidas em relacao a [HA2003].
            % NAO e um teste de igualdade estrita com o artigo.
            caso = main_hadi_nao_linear('caso');
            [w, S, u] = fem_nao_linear_solver(caso, caso.ref_areas);

            max_u = max(abs(u));
            max_S = max(abs(S));

            fprintf(['\n[CARACTERIZACAO HADI Case 2] solucao de referencia do artigo\n' ...
                     '  peso calculado : %9.4f kg   (artigo: %.1f kg, dif %+.3f%%)\n' ...
                     '  max |u|        : %9.4f mm   (limite: %.2f mm)\n' ...
                     '  max |sigma|    : %9.4f MPa  (limite: %.2f MPa)\n'], ...
                    w, caso.ref_peso, 100*(w-caso.ref_peso)/caso.ref_peso, ...
                    max_u, caso.d_max, max_S, caso.sigma_max);

            % Peso: proximo do publicado, mas nao identico (ver relatorio).
            tc.verifyEqual(w, caso.ref_peso, 'RelTol', 1e-3, ...
                'Peso afastou-se mais de 0.1% do valor publicado.');

            % A tensao respeita o limite do artigo.
            tc.verifyLessThan(max_S, caso.sigma_max, ...
                'Tensao maxima excedeu o admissivel.');

            % DISCREPANCIA DOCUMENTADA: o deslocamento excede levemente o
            % limite (51.02 mm > 50.80 mm). Estavel a variacoes de n_inc e
            % da magnitude exata da carga — nao e erro de convergencia.
            % Ver testHadi_UnidadesImperiaisExplicamOPeso abaixo.
            tc.verifyLessThan(max_u, 1.01 * caso.d_max, ...
                'Deslocamento excedeu o limite em mais de 1% (discrepancia mudou).');
        end

        function testHadi_UnidadesImperiaisExplicamOPeso(tc)
            % EXPLICACAO DA DISCREPANCIA DE PESO (investigacao de 2026-08-31).
            %
            % [HA2003] enuncia o problema em SI, mas os valores sao conversoes
            % de numeros redondos imperiais. A Tabela 1 foi calculada com
            % rho = 0.1 lb/in^3 EXATO (= 2767.990 kg/m^3), e nao com o
            % 2770 kg/m^3 declarado no texto.
            %
            % Este teste trava a explicacao: trocando SO a densidade pelo
            % valor imperial exato, o peso publicado e reproduzido.
            caso = main_hadi_nao_linear('caso');

            w_texto = fem_nao_linear_solver(caso, caso.ref_areas);

            c = caso;  c.dens = caso.dens_exata;
            w_exato = fem_nao_linear_solver(c, caso.ref_areas);

            fprintf(['\n[UNIDADES IMPERIAIS] peso da solucao de referencia\n' ...
                     '  dens = 2770.000 (texto do artigo) : %9.4f kg (%+.4f%%)\n' ...
                     '  dens = %7.3f (0.1 lb/in^3)      : %9.4f kg (%+.4f%%)\n'], ...
                    w_texto, 100*(w_texto-caso.ref_peso)/caso.ref_peso, ...
                    caso.dens_exata*1e9, ...
                    w_exato, 100*(w_exato-caso.ref_peso)/caso.ref_peso);

            % A densidade exata reproduz o publicado dentro do arredondamento
            % de impressao do proprio artigo (1 casa decimal em 2325.2).
            tc.verifyEqual(w_exato, caso.ref_peso, 'AbsTol', 0.1, ...
                'Densidade imperial exata deveria reproduzir a Tabela 1.');

            % E melhora em mais de uma ordem de grandeza sobre o valor do texto.
            err_texto = abs(w_texto - caso.ref_peso);
            err_exato = abs(w_exato - caso.ref_peso);
            tc.verifyLessThan(err_exato, err_texto/10, ...
                'A explicacao imperial deveria reduzir o erro em 10x ou mais.');
        end

        function testHadi_TodasAsSolucoesPublicadasReproduzemOPeso(tc)
            % A explicacao imperial nao pode valer so para uma linha da
            % Tabela 1. Este teste a aplica a TODAS as colunas do artigo que
            % listam areas — inclusive as de Gutkowski-Zawidzka [4], obtidas
            % por outro metodo. Se o peso de todas bate, a formula do peso e
            % o conjunto de comprimentos estao certos.
            caso = main_hadi_nao_linear('caso');
            c = caso;  c.dens = caso.dens_exata;

            % [HA2003] Tabela 1 — {rotulo, areas [mm^2], peso publicado [kg]}
            T = { ...
              'C1 G-Z Enumeration  ', [23226 323 17419 12258 323 1290 4516 12258 12258 65], 2429.5; ...
              'C1 Hadi discreto    ', [17419 65 17419 12258 65 65 4516 12258 17419 65],     2423.2; ...
              'C2 G-Z Sequential   ', [19355 65 19355 9677 65 65 5161 12903 12903 65],      2340.4; ...
              'C2 G-Z Enumeration  ', [19355 65 16129 7742 65 645 5161 12903 16129 65],     2339.9; ...
              'C2 Hadi discreto    ', caso.ref_areas,                                       2325.2 };

            fprintf('\n[TABELA 1 COMPLETA] peso com densidade imperial exata\n');
            for k = 1:size(T,1)
                w = fem_nao_linear_solver(c, T{k,2});
                fprintf('  %s: %9.2f kg (publicado %7.1f, %+.4f%%)\n', ...
                        T{k,1}, w, T{k,3}, 100*(w-T{k,3})/T{k,3});
                tc.verifyEqual(w, T{k,3}, 'RelTol', 1e-4, ...
                    sprintf('Peso de "%s" nao reproduz a Tabela 1.', strtrim(T{k,1})));
            end
        end

        function testHadi_NumeracaoDasDiagonaisEhAMelhorPossivel(tc)
            % As barras 7..10 sao as diagonais e tem TODAS o mesmo
            % comprimento (L*sqrt(2)). Logo permutar areas entre elas NAO
            % altera o peso — so o deslocamento. Isso levantou a hipotese de
            % que a numeracao das diagonais em [HA2003] Fig. 1 fosse outra, o
            % que explicaria o deslocamento acima do limite.
            %
            % HIPOTESE REFUTADA, e este teste trava a refutacao: a ordem
            % adotada aqui produz o MENOR deslocamento entre as 24
            % permutacoes possiveis. Qualquer outra numeracao seria pior, e
            % nenhuma torna a solucao de referencia viavel.
            caso = main_hadi_nao_linear('caso');
            c = caso;  c.E = caso.E_exato;  c.dens = caso.dens_exata;
            ref = caso.ref_areas;

            % Confirma a premissa: as quatro diagonais tem comprimento igual
            L = zeros(1,10);
            for e = 1:10
                n1 = caso.elements(1,e);  n2 = caso.elements(2,e);
                L(e) = norm(caso.nodes0(:,n2) - caso.nodes0(:,n1));
            end
            tc.verifyEqual(L(8:10), L(7)*ones(1,3), 'RelTol', 1e-12, ...
                'As barras 7..10 deveriam ter comprimento igual.');

            [~,~,u_ref] = fem_nao_linear_solver(c, ref);
            d_ref = max(abs(u_ref));

            P = perms(7:10);
            d_min_alt = inf;
            for k = 1:size(P,1)
                a = ref;  a(7:10) = ref(P(k,:));
                if isequal(a, ref), continue; end
                [w_p,~,u_p] = fem_nao_linear_solver(c, a);
                % O peso e invariante a permutacao das diagonais
                tc.verifyEqual(w_p, fem_nao_linear_solver(c, ref), 'RelTol', 1e-12);
                d_min_alt = min(d_min_alt, max(abs(u_p)));
            end

            fprintf(['\n[NUMERACAO DAS DIAGONAIS] max|u| da solucao de referencia\n' ...
                     '  ordem adotada        : %8.3f mm\n' ...
                     '  melhor alternativa   : %8.3f mm\n' ...
                     '  limite [HA2003]      : %8.3f mm\n'], ...
                    d_ref, d_min_alt, caso.d_max);

            tc.verifyLessThan(d_ref, d_min_alt, ...
                ['A numeracao adotada deixou de ser a melhor das 24 ' ...
                 'permutacoes das diagonais — a topologia mudou.']);
            tc.verifyGreaterThan(d_min_alt, caso.d_max, ...
                'Alguma permutacao tornou a referencia viavel — reabrir a hipotese.');
        end

        function testPesoNaoDependeDoCarregamento(tc)
            % O peso e propriedade do projeto (area x comprimento inicial x
            % densidade) e nao pode variar com a carga aplicada.
            caso = main_hadi_nao_linear('caso');

            w1 = fem_nao_linear_solver(caso, caso.ref_areas);

            prob = caso;
            prob.F_total = caso.F_total * 0.3;
            w2 = fem_nao_linear_solver(prob, caso.ref_areas);

            tc.verifyEqual(w1, w2, 'RelTol', 1e-12, ...
                'O peso nao pode depender do carregamento.');
        end

        function testPesoConfereComCalculoManual(tc)
            % Verificacao independente: soma de A_i * L0_i * densidade.
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            n_el = size(caso.elements, 2);
            w_manual = 0;
            for n = 1:n_el
                ni = caso.elements(1,n);  nj = caso.elements(2,n);
                dx = caso.nodes0(1,nj) - caso.nodes0(1,ni);
                dy = caso.nodes0(2,nj) - caso.nodes0(2,ni);
                w_manual = w_manual + areas(n) * sqrt(dx^2+dy^2) * caso.dens;
            end

            w_solver = fem_nao_linear_solver(caso, areas);
            tc.verifyEqual(w_solver, w_manual, 'RelTol', 1e-12);
        end

        % -----------------------------------------------------------------
        % 6. VALIDACAO DE ENTRADA
        % -----------------------------------------------------------------
        function testRejeitaVetorDeAreasComTamanhoErrado(tc)
            caso = main_hadi_nao_linear('caso');
            tc.verifyError(@() fem_nao_linear_solver(caso, ones(1,9)), ...
                'fem_nao_linear_solver:areasIncompativeis');
        end

        function testRejeitaAreaNaoPositiva(tc)
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;  areas(3) = 0;
            tc.verifyError(@() fem_nao_linear_solver(caso, areas), ...
                'fem_nao_linear_solver:areaNaoPositiva');
        end

        function testRejeitaProblemaIncompleto(tc)
            caso = main_hadi_nao_linear('caso');
            caso = rmfield(caso, 'E');
            tc.verifyError(@() fem_nao_linear_solver(caso, ones(1,10)), ...
                'fem_nao_linear_solver:campoAusente');
        end

        % -----------------------------------------------------------------
        % 7. EQUIVALENCIA COM O SOLVER ORIGINAL (arquivado em 05_legado)
        % -----------------------------------------------------------------
        function testSolverOriginalProduzResultadoIdentico(tc)
            % O solver ORIGINAL (pre-Bloco 2), arquivado em 05_legado/,
            % deve devolver exatamente o mesmo que o solver generico,
            % garantindo que a migracao nao alterou nenhum resultado.
            caso  = main_hadi_nao_linear('caso');
            areas = caso.ref_areas;

            [w1, S1, u1] = fem_truss_nonlinear_pre_bloco2(areas);
            [w2, S2, u2] = fem_nao_linear_solver(caso, areas);

            tc.verifyEqual(w1, w2, 'AbsTol', 0, 'Peso divergiu do solver original.');
            tc.verifyEqual(S1, S2, 'AbsTol', 0, 'Tensoes divergiram do solver original.');
            tc.verifyEqual(u1, u2, 'AbsTol', 0, 'Deslocamentos divergiram do solver original.');
        end
    end

    % ---------------------------------------------------------------------
    methods (Access = private)
        function F_int = montarForcasInternas(~, caso, d)
            % Reconstroi o vetor de forcas internas a partir das forcas
            % axiais e da geometria deformada — calculo independente do que
            % o solver faz internamente.
            s_dof = 2 * size(caso.nodes0, 2);
            F_int = zeros(s_dof, 1);
            for n = 1:size(caso.elements, 2)
                ni = caso.elements(1,n);  nj = caso.elements(2,n);
                dxc = d.nodes_cur(1,nj) - d.nodes_cur(1,ni);
                dyc = d.nodes_cur(2,nj) - d.nodes_cur(2,ni);
                Ln  = sqrt(dxc^2 + dyc^2);
                lam = dxc/Ln;  mu = dyc/Ln;
                idx = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
                F_int(idx) = F_int(idx) + d.P_axial(n) * [-lam; -mu; lam; mu];
            end
        end
    end
end
