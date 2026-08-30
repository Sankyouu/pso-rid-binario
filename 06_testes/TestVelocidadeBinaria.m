classdef TestVelocidadeBinaria < matlab.unittest.TestCase
% TESTVELOCIDADEBINARIA  Testes de rid_velocidade_binaria contra as Tabelas 1 e 2.
%
% BLOCO TESTADO: 1 (PSO-RID). Puramente logico — nenhum FEM envolvido.
%
% REFERENCIA:
%   [DF2011] Datta, D.; Figueira, J.R. Applied Soft Computing 11 (2011)
%            3625-3633. Secao 4.3, Eq. (7), Eq. (8), Tabelas 1 e 2.
%
% Eq. (7): transicoes admissiveis
%     x = 0  =>  v pertence a { 0, +1}      (Tabela 1)
%     x = 1  =>  v pertence a { 0, -1}      (Tabela 2)
%
% Eq. (8): mapeamento do valor real de velocidade
%     v > 0  =>  +1        v < 0  =>  -1
%     ambiguidade resolvida com probabilidade de mutacao p_m
%
% [D3][D4] Estes testes travam a correcao das duas divergencias em relacao
% a versao anterior do solver, que usava regra SIMETRICA v em {-1,0,+1}
% independentemente do bit atual.
%
% Executar:  runtests('06_testes/TestVelocidadeBinaria.m')

    properties (Constant)
        N_AMOSTRAS = 4000;   % amostras para os testes estatisticos
    end

    methods (TestClassSetup)
        function adicionarCaminhos(~)
            aqui = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(aqui, '..', '01_pso_rid')));
        end
    end

    methods (Test)

        % =================================================================
        % Eq. (7) — DOMINIO ADMISSIVEL DA VELOCIDADE
        % =================================================================
        function testEq7_ParaBitZeroVelocidadeNuncaEhNegativa(tc)
            % Tabela 1 [DF2011]: com x=0 os unicos valores admissiveis de v
            % sao 0 e +1. v = -1 e INADMISSIVEL (levaria a x = -1).
            %
            % [D3] REGRESSAO: a versao anterior podia gerar v = -1 aqui e
            % armazena-lo, corrompendo a inercia da iteracao seguinte.
            tendencias = [-1e6, -10, -1, -1e-9, 0, 1e-9, 1, 10, 1e6];
            for t = tendencias
                for pm = [0, 0.15, 0.5, 1.0]
                    for k = 1:40
                        v = rid_velocidade_binaria(0, t, pm);
                        tc.verifyTrue(v == 0 || v == 1, ...
                            sprintf(['Eq. (7)/Tab. 1 violada: x=0, tendencia=%g, ' ...
                                     'pm=%g gerou v=%g (esperado 0 ou +1).'], t, pm, v));
                    end
                end
            end
        end

        function testEq7_ParaBitUmVelocidadeNuncaEhPositiva(tc)
            % Tabela 2 [DF2011]: com x=1 os unicos valores admissiveis de v
            % sao 0 e -1. v = +1 e INADMISSIVEL (levaria a x = 2).
            tendencias = [-1e6, -10, -1, -1e-9, 0, 1e-9, 1, 10, 1e6];
            for t = tendencias
                for pm = [0, 0.15, 0.5, 1.0]
                    for k = 1:40
                        v = rid_velocidade_binaria(1, t, pm);
                        tc.verifyTrue(v == 0 || v == -1, ...
                            sprintf(['Eq. (7)/Tab. 2 violada: x=1, tendencia=%g, ' ...
                                     'pm=%g gerou v=%g (esperado 0 ou -1).'], t, pm, v));
                    end
                end
            end
        end

        function testEq7_PosicaoResultanteSempreBinaria(tc)
            % Consequencia direta da Eq. (7): x + v pertence sempre a {0,1},
            % sem necessidade de saturacao posterior.
            for x = [0 1]
                for t = [-1e3, -1, 0, 1, 1e3]
                    for pm = [0, 0.15, 1.0]
                        for k = 1:40
                            [~, x_next] = rid_velocidade_binaria(x, t, pm);
                            tc.verifyTrue(x_next == 0 || x_next == 1, ...
                                sprintf('Posicao resultante nao binaria: %g.', x_next));
                        end
                    end
                end
            end
        end

        function testCoerenciaEntreSaidas(tc)
            % x_next deve ser exatamente x_atual + v_next.
            for x = [0 1]
                for t = [-5, 0, 5]
                    for k = 1:100
                        [v, x_next] = rid_velocidade_binaria(x, t, 0.15);
                        tc.verifyEqual(x_next, x + v, ...
                            'x_next deve ser igual a x_atual + v_next.');
                    end
                end
            end
        end

        % =================================================================
        % Eq. (8) — MAPEAMENTO DETERMINISTICO (com pm = 0)
        % =================================================================
        function testEq8_SemMutacao_TendenciaPositivaComBitZeroLevaAUm(tc)
            % Tabela 1, tendencia > 0 => v = +1 (deterministico com pm=0).
            for t = [1e-9, 0.5, 1, 100, 1e6]
                for k = 1:50
                    [v, x_next] = rid_velocidade_binaria(0, t, 0);
                    tc.verifyEqual(v, 1, ...
                        sprintf('Eq. (8): x=0, tendencia=%g deveria dar v=+1.', t));
                    tc.verifyEqual(x_next, 1);
                end
            end
        end

        function testEq8_SemMutacao_TendenciaNegativaComBitZeroMantemZero(tc)
            % Tabela 1, tendencia < 0: v = -1 seria inadmissivel, logo v = 0.
            for t = [-1e-9, -0.5, -1, -100, -1e6]
                for k = 1:50
                    [v, x_next] = rid_velocidade_binaria(0, t, 0);
                    tc.verifyEqual(v, 0, ...
                        sprintf('Eq. (8): x=0, tendencia=%g deveria dar v=0.', t));
                    tc.verifyEqual(x_next, 0);
                end
            end
        end

        function testEq8_SemMutacao_TendenciaNegativaComBitUmLevaAZero(tc)
            % Tabela 2, tendencia < 0 => v = -1.
            for t = [-1e-9, -0.5, -1, -100, -1e6]
                for k = 1:50
                    [v, x_next] = rid_velocidade_binaria(1, t, 0);
                    tc.verifyEqual(v, -1, ...
                        sprintf('Eq. (8): x=1, tendencia=%g deveria dar v=-1.', t));
                    tc.verifyEqual(x_next, 0);
                end
            end
        end

        function testEq8_SemMutacao_TendenciaPositivaComBitUmMantemUm(tc)
            % Tabela 2, tendencia > 0: v = +1 seria inadmissivel, logo v = 0.
            for t = [1e-9, 0.5, 1, 100, 1e6]
                for k = 1:50
                    [v, x_next] = rid_velocidade_binaria(1, t, 0);
                    tc.verifyEqual(v, 0, ...
                        sprintf('Eq. (8): x=1, tendencia=%g deveria dar v=0.', t));
                    tc.verifyEqual(x_next, 1);
                end
            end
        end

        % =================================================================
        % COLUNA 7 DAS TABELAS 1 e 2 — PROBABILIDADE DE MUTACAO p_m
        % =================================================================
        function testPm_Zero_EhTotalmenteDeterministico(tc)
            % Com pm = 0 nunca se adota o valor alternativo.
            v_ref = rid_velocidade_binaria(0, 1.0, 0);
            for k = 1:200
                tc.verifyEqual(rid_velocidade_binaria(0, 1.0, 0), v_ref);
            end
        end

        function testPm_Um_SempreAdotaAlternativa(tc)
            % Com pm = 1 sempre se adota o valor ALTERNATIVO ao mapeado.
            % x=0, tendencia>0: mapeado = +1, logo alternativa = 0.
            for k = 1:200
                tc.verifyEqual(rid_velocidade_binaria(0, 1.0, 1.0), 0, ...
                    'Com pm=1 deveria adotar sempre a alternativa.');
            end
            % x=1, tendencia<0: mapeado = -1, logo alternativa = 0.
            for k = 1:200
                tc.verifyEqual(rid_velocidade_binaria(1, -1.0, 1.0), 0);
            end
        end

        function testPm_FrequenciaDaAlternativaBateComOValorNominal(tc)
            % A fracao de vezes em que a alternativa e adotada deve convergir
            % para p_m. Teste estatistico com tolerancia generosa.
            rng(20260828);   % reprodutibilidade
            for pm = [0.05, 0.15, 0.30]
                n_alt = 0;
                for k = 1:tc.N_AMOSTRAS
                    v = rid_velocidade_binaria(0, 1.0, pm);   % mapeado = +1
                    if v == 0, n_alt = n_alt + 1; end          % alternativa = 0
                end
                freq = n_alt / tc.N_AMOSTRAS;
                tc.verifyEqual(freq, pm, 'AbsTol', 0.03, ...
                    sprintf('Frequencia da alternativa (%.3f) longe de pm=%.3f.', freq, pm));
            end
        end

        % =================================================================
        % EMPATE EXATO (tendencia == 0)
        % =================================================================
        function testTendenciaZero_DivideMeioAMeio(tc)
            % Nas linhas das Tabelas 1 e 2 sem relacao definida, cada opcao
            % vale 50% ("In the case of double options, each value is
            % considered with 50% probability", [DF2011] Sec. 4.3).
            rng(20260828);
            n_um = 0;
            for k = 1:tc.N_AMOSTRAS
                if rid_velocidade_binaria(0, 0, 0) == 1
                    n_um = n_um + 1;
                end
            end
            freq = n_um / tc.N_AMOSTRAS;
            tc.verifyEqual(freq, 0.5, 'AbsTol', 0.03, ...
                sprintf('Empate deveria dividir 50/50, obteve %.3f.', freq));
        end

        % =================================================================
        % VALIDACAO DE ENTRADA
        % =================================================================
        function testRejeitaPosicaoNaoBinaria(tc)
            tc.verifyError(@() rid_velocidade_binaria(2,   1, 0.15), ...
                'rid_velocidade_binaria:posicaoInvalida');
            tc.verifyError(@() rid_velocidade_binaria(-1,  1, 0.15), ...
                'rid_velocidade_binaria:posicaoInvalida');
            tc.verifyError(@() rid_velocidade_binaria(0.5, 1, 0.15), ...
                'rid_velocidade_binaria:posicaoInvalida');
        end
    end
end
