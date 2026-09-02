classdef TestRegrasDeb < matlab.unittest.TestCase
% TESTREGRASDEB  Testes de domina_deb (regra de viabilidade de Deb).
%
% BLOCO TESTADO: 1 (PSO-RID). Puramente logico — nenhum FEM envolvido.
%
% REFERENCIA:
%   [DEB2000] Deb, K. "An efficient constraint handling method for genetic
%             algorithms". Comput. Methods Appl. Mech. Engrg. 186 (2000)
%             311-338.  ->  pag. 316
%
%   Citacao literal dos tres criterios:
%     "1. Any feasible solution is preferred to any infeasible solution.
%      2. Among two feasible solutions, the one having better objective
%         function value is preferred.
%      3. Among two infeasible solutions, the one having smaller constraint
%         violation is preferred."
%
% Executar:  runtests('06_testes/TestRegrasDeb.m')

    properties
        % Handles das funcoes locais de pso_rid.m (ver "ACESSO AOS
        % AUXILIARES PARA TESTE" no cabecalho daquele arquivo).
        aux
    end

    methods (TestClassSetup)
        function adicionarCaminhos(tc)
            aqui = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(aqui, '..', '01_pso_rid')));
            tc.aux = pso_rid('auxiliares');
        end
    end

    methods (Test)

        % -----------------------------------------------------------------
        % CRITERIO 1 — "Any feasible solution is preferred to any infeasible"
        % -----------------------------------------------------------------
        function testCriterio1_ViavelDominaInviavel(tc)
            % Viavel com custo PESSIMO ainda domina inviavel com custo otimo.
            % Esse e o ponto central da abordagem sem parametro de penalizacao.
            tc.verifyTrue(tc.aux.domina_deb(1e9, 0, 1.0, 5.0), ...
                'Criterio 1: viavel (custo alto) deveria dominar inviavel.');
        end

        function testCriterio1_InviavelNaoDominaViavel(tc)
            tc.verifyFalse(tc.aux.domina_deb(1.0, 5.0, 1e9, 0), ...
                'Criterio 1: inviavel jamais deveria dominar viavel.');
        end

        function testCriterio1_ViolacaoMinusculaAindaEhInviavel(tc)
            % Qualquer violacao > 0 torna a solucao inviavel.
            tc.verifyTrue(tc.aux.domina_deb(1e6, 0, 1.0, 1e-12), ...
                'Violacao positiva minuscula ainda caracteriza inviabilidade.');
        end

        % -----------------------------------------------------------------
        % CRITERIO 2 — "Among two feasible, better objective is preferred"
        % -----------------------------------------------------------------
        function testCriterio2_EntreViaveisVenceMenorCusto(tc)
            tc.verifyTrue (tc.aux.domina_deb(10, 0, 20, 0));
            tc.verifyFalse(tc.aux.domina_deb(20, 0, 10, 0));
        end

        function testCriterio2_EmpateDeCustoNaoDomina(tc)
            % Estritamente melhor: empate nao substitui a referencia
            % (evita substituicoes desnecessarias e churn no p-best/g-best).
            tc.verifyFalse(tc.aux.domina_deb(10, 0, 10, 0), ...
                'Empate de custo entre viaveis nao deveria dominar.');
        end

        % -----------------------------------------------------------------
        % CRITERIO 3 — "Among two infeasible, smaller violation is preferred"
        % -----------------------------------------------------------------
        function testCriterio3_EntreInviaveisVenceMenorViolacao(tc)
            tc.verifyTrue (tc.aux.domina_deb(1e9, 2.0, 1.0, 5.0), ...
                'Criterio 3: menor violacao vence, INDEPENDENTE do custo.');
            tc.verifyFalse(tc.aux.domina_deb(1.0, 5.0, 1e9, 2.0));
        end

        function testCriterio3_CustoEhIrrelevanteEntreInviaveis(tc)
            % [DEB2000] pag. 316: entre inviaveis compara-se SOMENTE a
            % violacao. O custo nao deve influenciar em nada.
            tc.verifyTrue(tc.aux.domina_deb(0,    3.0, 1e300, 4.0));
            tc.verifyTrue(tc.aux.domina_deb(1e300, 3.0, 0,     4.0), ...
                'Entre inviaveis, o custo nao pode influenciar a decisao.');
        end

        function testCriterio3_EmpateDeViolacaoNaoDomina(tc)
            tc.verifyFalse(tc.aux.domina_deb(1, 5.0, 1e9, 5.0), ...
                'Empate de violacao entre inviaveis nao deveria dominar.');
        end

        % -----------------------------------------------------------------
        % INICIALIZACAO (referencia = Inf/Inf, como no inicio do PSO)
        % -----------------------------------------------------------------
        function testQualquerCoisaDominaReferenciaInicial(tc)
            % No inicio, p-best e g-best comecam com custo=Inf e violacao=Inf.
            % Qualquer candidato — viavel ou nao — deve substitui-los.
            tc.verifyTrue(tc.aux.domina_deb(100, 0,   inf, inf), ...
                'Candidato viavel deveria dominar a referencia inicial.');
            tc.verifyTrue(tc.aux.domina_deb(100, 50,  inf, inf), ...
                'Candidato inviavel deveria dominar a referencia inicial.');
        end

        function testReferenciaInicialNaoDominaNada(tc)
            tc.verifyFalse(tc.aux.domina_deb(inf, inf, 100, 0));
            tc.verifyFalse(tc.aux.domina_deb(inf, inf, 100, 50));
        end

        % -----------------------------------------------------------------
        % ESTOURO DE CATALOGO (violacao = Inf)
        % -----------------------------------------------------------------
        function testQualquerProjetoAvaliavelVenceUmEstouro(tc)
            % Uma particula fora do catalogo nao e um projeto ruim: e a
            % AUSENCIA de projeto. O solver lhe atribui violacao = Inf, e dai
            % as tres regras garantem que ela perca para qualquer projeto
            % avaliavel, por pior que ele seja.
            %
            % Por que Inf e nao uma constante grande: violacoes reais medidas
            % neste problema chegam a 4.5e8, entao qualquer constante fixa
            % seria ultrapassada e o criterio 3 voltaria a preferir o estouro.
            for viol_real = [2e-5, 1, 9.1e4, 4.5e8, 1e300]
                tc.verifyTrue(tc.aux.domina_deb(2500, viol_real, inf, inf), ...
                    sprintf(['Projeto avaliavel com violacao %g deveria vencer ' ...
                             'um estouro de catalogo.'], viol_real));
                tc.verifyFalse(tc.aux.domina_deb(inf, inf, 2500, viol_real), ...
                    'Estouro nunca deveria vencer um projeto avaliavel.');
            end
        end

        function testNaoHaOrdenacaoEntreEstouros(tc)
            % Consequencia deliberada de usar Inf: todo estouro colapsa no
            % mesmo ponto, entao nenhum domina outro. A pressao que traz a
            % particula de volta ao catalogo vem do g-best pela Eq. (1), nao
            % da comparacao entre posicoes sem significado fisico.
            tc.verifyFalse(tc.aux.domina_deb(inf, inf, inf, inf), ...
                'Um estouro nao deveria dominar outro estouro.');
        end

        % -----------------------------------------------------------------
        % CONVENCAO DE VIABILIDADE
        % -----------------------------------------------------------------
        function testViolacaoZeroEhViavel(tc)
            tc.verifyTrue(tc.aux.domina_deb(1, 0, 2, 0), ...
                'Violacao exatamente 0 deve ser tratada como viavel.');
        end

        function testViolacaoNegativaEhViavel(tc)
            % Folga negativa (margem) tambem conta como viavel.
            tc.verifyTrue(tc.aux.domina_deb(1, -3, 2, 0));
            tc.verifyTrue(tc.aux.domina_deb(1, 0, 2, -3));
        end

        % -----------------------------------------------------------------
        % PROPRIEDADES ESTRUTURAIS DA RELACAO
        % -----------------------------------------------------------------
        function testAntissimetria(tc)
            % Para pares distintos, no maximo um dos sentidos pode dominar.
            pares = { 10,0, 20,0 ;      % ambos viaveis
                      10,0, 20,5 ;      % viavel vs inviavel
                      10,3, 20,5 ;      % ambos inviaveis
                      10,5, 20,3 };
            for i = 1:size(pares,1)
                a_c = pares{i,1}; a_v = pares{i,2};
                b_c = pares{i,3}; b_v = pares{i,4};
                ab = tc.aux.domina_deb(a_c, a_v, b_c, b_v);
                ba = tc.aux.domina_deb(b_c, b_v, a_c, a_v);
                tc.verifyFalse(ab && ba, ...
                    sprintf('Antissimetria violada no par %d.', i));
            end
        end

        function testIrreflexividade(tc)
            % Nenhuma solucao domina a si mesma.
            tc.verifyFalse(tc.aux.domina_deb(10, 0,  10, 0));
            tc.verifyFalse(tc.aux.domina_deb(10, 5,  10, 5));
        end

        function testTransitividadeEmAmostraExaustiva(tc)
            % Varre uma grade de (custo, violacao) e confere transitividade:
            % se A domina B e B domina C, entao A domina C.
            custos    = [1 5 10];
            violacoes = [0 0 2 7];
            solucoes  = [];
            for c = custos
                for v = violacoes
                    solucoes(end+1,:) = [c v]; %#ok<AGROW>
                end
            end

            n = size(solucoes,1);
            for a = 1:n
                for b = 1:n
                    if ~tc.aux.domina_deb(solucoes(a,1),solucoes(a,2), ...
                                       solucoes(b,1),solucoes(b,2))
                        continue;
                    end
                    for c = 1:n
                        if tc.aux.domina_deb(solucoes(b,1),solucoes(b,2), ...
                                          solucoes(c,1),solucoes(c,2))
                            tc.verifyTrue( ...
                                tc.aux.domina_deb(solucoes(a,1),solucoes(a,2), ...
                                               solucoes(c,1),solucoes(c,2)), ...
                                sprintf(['Transitividade violada: ' ...
                                         '(%g,%g) > (%g,%g) > (%g,%g).'], ...
                                        solucoes(a,:), solucoes(b,:), solucoes(c,:)));
                        end
                    end
                end
            end
        end
    end
end
