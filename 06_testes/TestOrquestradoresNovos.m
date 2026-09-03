classdef TestOrquestradoresNovos < matlab.unittest.TestCase
% TESTORQUESTRADORESNOVOS  Reproducao dos valores publicados nos 4 casos novos.
%
% CASOS COBERTOS
%   main_hadi_20barras      [HA2003] Sec. 2.2  — trelica de 20 barras
%   main_hadi_51barras      [HA2003] Sec. 6.3  — trelica de teto, 51 barras
%   main_datta_engrenagens  [DF2011] Sec. 5.1  — trem de engrenagens
%   main_datta_mola         [DF2011] Sec. 5.2  — mola helicoidal
%
% O QUE ESTA CLASSE PROTEGE
%
% Cada um dos quatro foi conferido UMA vez, na mao, contra os valores
% publicados. Sem teste, uma regressao no pso_rid ou no solver FEM quebraria
% essa concordancia sem que nada acusasse. Aqui a concordancia fica travada.
%
% Nao se testa a QUALIDADE da otimizacao — isso e papel do estudo
% estatistico. Testa-se que, dada a solucao publicada, este codigo devolve o
% valor publicado. E um teste de modelagem, nao de convergencia, e por isso
% roda em segundos: nenhuma execucao do PSO e disparada.
%
% Executar:  runtests('06_testes/TestOrquestradoresNovos.m')

    properties
        aux20   % auxiliares de main_hadi_20barras
        aux51   % auxiliares de main_hadi_51barras
        auxEng  % auxiliares de main_datta_engrenagens
        auxMola % auxiliares de main_datta_mola
    end

    methods (TestClassSetup)
        function adicionarCaminhos(tc)
            aqui = fileparts(mfilename('fullpath'));
            raiz = fullfile(aqui, '..');
            addpath(raiz);
            addpath(genpath(fullfile(raiz, '01_pso_rid')));
            addpath(genpath(fullfile(raiz, '02_fem_nao_linear')));
            addpath(genpath(fullfile(raiz, '03_orquestrador')));

            tc.aux20   = main_hadi_20barras('auxiliares');
            tc.aux51   = main_hadi_51barras('auxiliares');
            tc.auxEng  = main_datta_engrenagens('auxiliares');
            tc.auxMola = main_datta_mola('auxiliares');
        end
    end

    methods (Test)

        % -----------------------------------------------------------------
        % 1. TRELICA DE 20 BARRAS — [HA2003] Sec. 2.2, Tabela 2
        % -----------------------------------------------------------------
        function test20Barras_ReproduzOsDoisPesosPublicados(tc)
            % Tabela 2 traz duas solucoes: Case 1 (nao linear, 6222 kg) e
            % Case 2 (linear, 6253 kg). Ambas precisam sair deste codigo.
            %
            % TOLERANCIA de 0.1%: a tabela publica os pesos arredondados ao
            % quilograma, o que ja vale ate 0.016% de incerteza em 6222 kg.

            caso = main_hadi_20barras('caso');

            peso_nl = tc.aux20.avaliar_projeto(caso.ref_areas_nao_linear, caso, 'nao_linear');
            peso_li = tc.aux20.avaliar_projeto(caso.ref_areas_linear,     caso, 'linear');

            tc.verifyEqual(peso_nl, caso.ref_peso_nao_linear, 'RelTol', 1e-3, ...
                'Case 1 divergiu do peso publicado (6222 kg).');
            tc.verifyEqual(peso_li, caso.ref_peso_linear, 'RelTol', 1e-3, ...
                'Case 2 divergiu do peso publicado (6253 kg).');

            fprintf(['\n[20 BARRAS] Case 1 %.1f kg (publicado %d) | ' ...
                     'Case 2 %.1f kg (publicado %d)\n'], ...
                    peso_nl, caso.ref_peso_nao_linear, ...
                    peso_li, caso.ref_peso_linear);
        end

        function test20Barras_AgrupamentoCobreTodasAsBarras(tc)
            % O PSO otimiza 10 grupos; o FEM precisa de 20 areas. A ponte e
            % expandir_areas. Se o mapa de grupos tiver furo, o FEM recebe
            % area zero e o resultado fica silenciosamente errado.
            caso = main_hadi_20barras('caso');

            areas_grupos = 1:numel(caso.ref_areas_nao_linear);
            areas_barras = tc.aux20.expandir_areas(areas_grupos, caso);

            tc.verifyEqual(numel(areas_barras), size(caso.elements, 2), ...
                'expandir_areas deve devolver uma area por barra.');
            tc.verifyTrue(all(areas_barras > 0), ...
                'Nenhuma barra pode ficar sem grupo.');
            tc.verifyEqual(unique(areas_barras(:))', areas_grupos, ...
                'Todo grupo deve ser usado por pelo menos uma barra.');
        end

        % -----------------------------------------------------------------
        % 2. TRELICA DE 51 BARRAS — [HA2003] Sec. 6.3, Tabela 3
        % -----------------------------------------------------------------
        function test51Barras_ReproduzVolumesDosCasos2e3(tc)
            % O objetivo aqui e VOLUME, nao peso (caso.dens = 1). Case 2 sem
            % flambagem, Case 3 com. As duas solucoes estao no catalogo, ou
            % seja, sao alvos legitimos deste codigo discreto.

            caso = main_hadi_51barras('caso');

            vol2 = tc.aux51.avaliar_projeto(caso.ref_areas_case2, caso, 'nao_linear', false);
            vol3 = tc.aux51.avaliar_projeto(caso.ref_areas_case3, caso, 'nao_linear', true);

            tc.verifyEqual(vol2, caso.ref_volume_case2, 'RelTol', 5e-4, ...
                'Case 2 divergiu do volume publicado.');
            tc.verifyEqual(vol3, caso.ref_volume_case3, 'RelTol', 5e-4, ...
                'Case 3 divergiu do volume publicado.');

            fprintf(['\n[51 BARRAS] Case 2 %.4e mm3 (publicado %.4e) | ' ...
                     'Case 3 %.4e mm3 (publicado %.4e)\n'], ...
                    vol2, caso.ref_volume_case2, vol3, caso.ref_volume_case3);
        end

        function test51Barras_DuasHipotesesDeCargaSaoResolvidas(tc)
            % Este e o unico caso do projeto com mais de uma hipotese de
            % carga. Se so a primeira for resolvida, a restricao da segunda
            % desaparece e o otimizador encontra pesos irreais.
            caso = main_hadi_51barras('caso');
            areas_barras = tc.aux51.expandir_areas(caso.ref_areas_case2, caso);

            [Sigma, u, L0] = tc.aux51.resolver_todas_hipoteses(areas_barras, caso, 'nao_linear');

            tc.verifyEqual(numel(Sigma), 2, ...
                'As duas hipoteses de carga precisam ser resolvidas.');
            tc.verifyEqual(numel(u), 2);
            tc.verifyEqual(numel(L0), size(caso.elements, 2));
            tc.verifyNotEqual(Sigma{1}, Sigma{2}, ...
                'Hipoteses distintas nao podem dar tensoes identicas.');
        end

        % -----------------------------------------------------------------
        % 3. TREM DE ENGRENAGENS — [DF2011] Sec. 5.1
        % -----------------------------------------------------------------
        function testEngrenagens_OtimoGlobalDaOValorPublicado(tc)
            % f* = 2.7009e-12 para z = (16,19,43,49). Aqui f_otimo e
            % calculado a partir de z_otimo, entao o teste verifica que
            % avaliar_projeto concorda com essa definicao — ou seja, que a
            % funcao objetivo usada pelo PSO e mesmo a Eq. do artigo.
            caso = main_datta_engrenagens('caso');

            [f, viol] = tc.auxEng.avaliar_projeto(caso.z_otimo, caso);

            tc.verifyEqual(f, caso.f_otimo, 'RelTol', 1e-12, ...
                'A funcao objetivo divergiu do otimo global tabelado.');
            tc.verifyEqual(viol, 0, ...
                'O otimo global e viavel: violacao deve ser zero.');
            tc.verifyEqual(tc.auxEng.razao(caso.z_otimo), (16*19)/(43*49), ...
                'RelTol', 1e-12, 'A razao de engrenagens saiu errada.');

            fprintf('\n[ENGRENAGENS] f(16,19,43,49) = %.4e (publicado 2.7009e-12)\n', f);
        end

        function testEngrenagens_SolucoesDaLiteraturaSaoPiores(tc)
            % As 6 solucoes publicadas por outros autores precisam dar f
            % MAIOR que o otimo global. Se alguma sair menor, a funcao
            % objetivo esta errada.
            caso = main_datta_engrenagens('caso');

            for i = 1:size(caso.ref_z, 1)
                fi = tc.auxEng.avaliar_projeto(caso.ref_z(i,:), caso);
                tc.verifyGreaterThanOrEqual(fi, caso.f_otimo, ...
                    sprintf('A solucao de %s ficou abaixo do otimo global.', ...
                            caso.ref_fontes{i}));
            end
        end

        % -----------------------------------------------------------------
        % 4. MOLA HELICOIDAL — [DF2011] Sec. 5.2
        % -----------------------------------------------------------------
        function testMola_OtimoPublicadoCaiNaFronteiraDeViabilidade(tc)
            % x = [D, N, d]. O D publicado tem 6 casas decimais, e a restricao
            % g8 esta ATIVA no otimo: arredondar D para baixo joga o ponto
            % 3.7e-08 para FORA da regiao viavel.
            %
            %   D publicado ....... 1.223041      -> violacao 2.44e-08
            %   fronteira real .... 1.223041037   (bissecao, 2026-09-02)
            %
            % Por isso nao se exige viabilidade estrita: exige-se que a
            % violacao seja da ordem do arredondamento publicado (< 1e-6) e
            % que o topo da faixa publicada ja seja estritamente viavel. Um
            % erro de modelagem de verdade sairia varias ordens acima disso.
            caso = main_datta_mola('caso');

            [f, viol, g] = tc.auxMola.avaliar_projeto( ...
                [caso.ref_D, caso.ref_N, caso.ref_d], caso);

            tc.verifyGreaterThanOrEqual(f, caso.ref_faixa_f(1), ...
                'f ficou abaixo da faixa publicada.');
            tc.verifyLessThanOrEqual(f, caso.ref_faixa_f(2), ...
                'f ficou acima da faixa publicada.');
            tc.verifyLessThan(viol, 1e-6, ...
                'A violacao passou da escala do arredondamento publicado.');
            tc.verifyEqual(find(g == max(g), 1), 8, ...
                'A restricao ativa no otimo publicado deve ser g8.');

            % O topo da faixa publicada de D e estritamente viavel.
            [~, viol_topo] = tc.auxMola.avaliar_projeto( ...
                [caso.ref_faixa_D(2), caso.ref_N, caso.ref_d], caso);
            tc.verifyEqual(viol_topo, 0, ...
                'O topo de ref_faixa_D deveria ser estritamente viavel.');

            fprintf('\n[MOLA] f = %.6f (publicado %.6f), violacao %.3e (g8 ativa)\n', ...
                    f, caso.ref_f, viol);
        end

        function testMola_SolucoesDaLiteraturaReproduzemOsFPublicados(tc)
            % caso.ref_x e caso.ref_f_pub vem da Tabela 6. Reproduzir os
            % sete valores de f a partir dos sete x e o teste mais forte de
            % que a funcao objetivo esta correta.
            caso = main_datta_mola('caso');

            for i = 1:size(caso.ref_x, 1)
                fi = tc.auxMola.avaliar_projeto(caso.ref_x(i,:), caso);
                tc.verifyEqual(fi, caso.ref_f_pub(i), 'RelTol', 1e-3, ...
                    sprintf('f divergiu do publicado para %s.', caso.ref_fontes{i}));
            end
        end

        % -----------------------------------------------------------------
        % 5. INTERFACE COM O OTIMIZADOR
        % -----------------------------------------------------------------
        function testTodosOsCasos_ExpoemConfigVarsValidoParaOPso(tc)
            % Qualquer caso precisa entregar config_vars que o pso_rid saiba
            % mapear. mapear_dimensoes falha alto se algo estiver errado.
            aux = pso_rid('auxiliares');

            casos = {main_hadi_20barras('caso'),     'main_hadi_20barras'; ...
                     main_hadi_51barras('caso'),     'main_hadi_51barras'; ...
                     main_datta_engrenagens('caso'), 'main_datta_engrenagens'; ...
                     main_datta_mola('caso'),        'main_datta_mola'};

            for k = 1:size(casos, 1)
                caso = casos{k,1};
                nome = casos{k,2};

                tc.verifyTrue(isfield(caso, 'config_vars'), ...
                    sprintf('%s nao expoe config_vars.', nome));

                for i = 1:numel(caso.config_vars)
                    tc.verifyTrue(any(strcmp(caso.config_vars(i).tipo, {'R','I','D'})), ...
                        sprintf('%s: tipo de variavel invalido em config_vars(%d).', nome, i));
                end

                [mapa, total_dim] = aux.mapear_dimensoes(caso.config_vars);
                tc.verifyEqual(numel(mapa), numel(caso.config_vars), ...
                    sprintf('%s: o mapa de dimensoes nao cobre todas as variaveis.', nome));
                tc.verifyGreaterThan(total_dim, 0, ...
                    sprintf('%s: dimensao total do enxame nao pode ser zero.', nome));
            end
        end

    end
end
