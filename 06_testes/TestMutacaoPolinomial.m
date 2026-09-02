classdef TestMutacaoPolinomial < matlab.unittest.TestCase
% TESTMUTACAOPOLINOMIAL  Testes de mutacao_polinomial (Eq. 9).
%
% BLOCO TESTADO: 1 (PSO-RID). Puramente numerico — nenhum FEM envolvido.
%
% REFERENCIA:
%   [DF2011] Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633,
%            Eq. (9), Secao 5.
%
%   x  <-  x + (x^u - x^l) * d_q
%
%   d_q = [2r + (1-2r)*((x^u - x)/(x^u - x^l))^(eta+1)]^(1/(eta+1)) - 1,  r < 0.5
%   d_q = 1 - [2(1-r) + (2r-1)*((x - x^l)/(x^u - x^l))^(eta+1)]^(1/(eta+1)), senao
%
% [D5] Substitui a perturbacao gaussiana usada na versao anterior.
%
% Executar:  runtests('06_testes/TestMutacaoPolinomial.m')

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
        % LIMITES DO DOMINIO
        % -----------------------------------------------------------------
        function testSempreDentroDosLimites(tc)
            rng(20260828);
            x_l = -3.5;  x_u = 7.25;
            for eta = [25, 35, 45]
                for k = 1:2000
                    x      = x_l + rand * (x_u - x_l);
                    x_novo = tc.aux.mutacao_polinomial(x, x_l, x_u, eta);
                    tc.verifyGreaterThanOrEqual(x_novo, x_l, ...
                        'Mutacao saiu abaixo do limite inferior.');
                    tc.verifyLessThanOrEqual(x_novo, x_u, ...
                        'Mutacao saiu acima do limite superior.');
                end
            end
        end

        function testNosExtremosDoDominioPermaneceValido(tc)
            rng(20260828);
            x_l = 0;  x_u = 1;
            for k = 1:500
                tc.verifyGreaterThanOrEqual(tc.aux.mutacao_polinomial(x_l, x_l, x_u, 30), x_l);
                tc.verifyLessThanOrEqual   (tc.aux.mutacao_polinomial(x_u, x_l, x_u, 30), x_u);
            end
        end

        % -----------------------------------------------------------------
        % PROPRIEDADE CENTRAL: eta maior => perturbacao menor
        % -----------------------------------------------------------------
        function testEtaMaiorProduzPerturbacaoMenor(tc)
            % Propriedade classica da mutacao polinomial: o indice de
            % distribuicao eta controla a concentracao em torno de x.
            % [DF2011] Sec. 5 sorteia eta em [25,45] — faixa de perturbacao
            % pequena, adequada a busca local (refinamento).
            rng(20260828);
            x_l = 0;  x_u = 100;  x = 50;
            n   = 4000;

            desvio = zeros(1,3);
            etas   = [5, 25, 45];
            for e = 1:numel(etas)
                acum = 0;
                for k = 1:n
                    acum = acum + abs(tc.aux.mutacao_polinomial(x, x_l, x_u, etas(e)) - x);
                end
                desvio(e) = acum / n;
            end

            tc.verifyGreaterThan(desvio(1), desvio(2), ...
                'eta=5 deveria perturbar mais que eta=25.');
            tc.verifyGreaterThan(desvio(2), desvio(3), ...
                'eta=25 deveria perturbar mais que eta=45.');

            fprintf(['\n[CARACTERIZACAO Eq.9] perturbacao media |x_novo - x| ' ...
                     'em dominio [0,100], x=50:\n' ...
                     '  eta= 5 -> %.3f\n  eta=25 -> %.3f\n  eta=45 -> %.3f\n'], ...
                    desvio(1), desvio(2), desvio(3));
        end

        function testPerturbacaoEhLocalNaFaixaDoArtigo(tc)
            % Com eta em [25,45] a mutacao deve ser de fato LOCAL: a
            % perturbacao media fica bem abaixo de 10% do dominio.
            rng(20260828);
            x_l = 0;  x_u = 100;  x = 50;  n = 4000;
            acum = 0;
            for k = 1:n
                eta  = 25 + rand * 20;             % faixa de [DF2011] Sec. 5
                acum = acum + abs(tc.aux.mutacao_polinomial(x, x_l, x_u, eta) - x);
            end
            media = acum / n;
            tc.verifyLessThan(media, 10, ...
                'Com eta em [25,45] a mutacao deveria ser local (<10% do dominio).');
        end

        % -----------------------------------------------------------------
        % SIMETRIA E CENTRAGEM
        % -----------------------------------------------------------------
        function testMutacaoEhCentradaNoValorAtual(tc)
            % Partindo do centro do dominio, a media das mutacoes deve ficar
            % proxima do proprio ponto (perturbacao aproximadamente simetrica).
            rng(20260828);
            x_l = 0;  x_u = 100;  x = 50;  n = 6000;
            acum = 0;
            for k = 1:n
                acum = acum + tc.aux.mutacao_polinomial(x, x_l, x_u, 30);
            end
            tc.verifyEqual(acum/n, 50, 'AbsTol', 1.5, ...
                'Mutacao a partir do centro deveria ser aproximadamente centrada.');
        end

        function testProduzValoresDosDoisLados(tc)
            % Deve gerar tanto valores acima quanto abaixo de x.
            rng(20260828);
            x_l = 0;  x_u = 10;  x = 5;
            acima = false;  abaixo = false;
            for k = 1:500
                x_novo = tc.aux.mutacao_polinomial(x, x_l, x_u, 30);
                if x_novo > x, acima  = true; end
                if x_novo < x, abaixo = true; end
            end
            tc.verifyTrue(acima,  'Nunca gerou valor acima de x.');
            tc.verifyTrue(abaixo, 'Nunca gerou valor abaixo de x.');
        end

        % -----------------------------------------------------------------
        % CASOS DEGENERADOS E VALIDACAO
        % -----------------------------------------------------------------
        function testDominioDegeneradoRetornaEntrada(tc)
            tc.verifyEqual(tc.aux.mutacao_polinomial(5, 5, 5, 30), 5, ...
                'Dominio de largura zero deveria devolver o proprio valor.');
            tc.verifyEqual(tc.aux.mutacao_polinomial(5, 7, 3, 30), 5, ...
                'Dominio invertido deveria devolver o proprio valor.');
        end

        function testRejeitaEtaNaoPositivo(tc)
            tc.verifyError(@() tc.aux.mutacao_polinomial(5, 0, 10, 0), ...
                'mutacao_polinomial:etaInvalido');
            tc.verifyError(@() tc.aux.mutacao_polinomial(5, 0, 10, -3), ...
                'mutacao_polinomial:etaInvalido');
        end

        function testDominioNegativoFunciona(tc)
            rng(20260828);
            x_l = -50;  x_u = -10;
            for k = 1:500
                x      = x_l + rand * (x_u - x_l);
                x_novo = tc.aux.mutacao_polinomial(x, x_l, x_u, 30);
                tc.verifyGreaterThanOrEqual(x_novo, x_l);
                tc.verifyLessThanOrEqual(x_novo, x_u);
            end
        end
    end
end
