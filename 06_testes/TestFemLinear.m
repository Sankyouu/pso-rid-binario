classdef TestFemLinear < matlab.unittest.TestCase
% TESTFEMLINEAR  Testes do solver elastico linear (Bloco 2a).
%
% BLOCO TESTADO: 2 (FEM linear), isoladamente.
%
% Estrategia: casos com solucao fechada calculavel a mao (barra unica sob
% tracao pura, trelica simetrica), equilibrio, simetria e propriedades
% estruturais da matriz de rigidez.
%
% Executar:  runtests('06_testes/TestFemLinear.m')

    methods (TestClassSetup)
        function adicionarCaminhos(~)
            aqui = fileparts(mfilename('fullpath'));
            raiz = fullfile(aqui, '..');
            addpath(genpath(fullfile(raiz, '01_pso_rid')));
            addpath(fullfile(raiz, '05_legado'));
            addpath(genpath(fullfile(raiz, '02_fem_nao_linear')));
        end
    end

    methods (Test)

        % =================================================================
        % 1. SOLUCAO FECHADA: BARRA UNICA SOB TRACAO PURA
        % =================================================================
        function testBarraUnica_DeslocamentoEhPLsobreEA(tc)
            % Caso mais elementar da resistencia dos materiais:
            %       delta = P*L / (E*A)      e      sigma = P/A
            % Calculavel a mao, sem nenhuma referencia externa.
            E = 200000;   % [MPa]
            A = 50;       % [mm^2]
            L = 2000;     % [mm]
            P = 10000;    % [N]

            prob.nodes0   = [0, L; 0, 0];    % no 1 na origem, no 2 a direita
            prob.elements = [1; 2];
            prob.apoios   = [1 2];           % no 1 totalmente fixo
            prob.E        = E;
            prob.dens     = 7.85e-6;
            prob.F_total  = [0; 0; P; 0];    % tracao no no 2, direcao x

            % O no 2 precisa de restricao vertical (barra horizontal nao tem
            % rigidez transversal): restringe o GDL 4.
            prob.apoios = [1 2 4];

            [~, Sigma, ~, d] = fem_linear_solver(prob, A);

            delta_esperado = P * L / (E * A);
            sigma_esperado = P / A;

            tc.verifyEqual(d.u(3), delta_esperado, 'RelTol', 1e-12, ...
                'Deslocamento axial diverge de P*L/(E*A).');
            tc.verifyEqual(Sigma(1), sigma_esperado, 'RelTol', 1e-12, ...
                'Tensao axial diverge de P/A.');
        end

        function testBarraUnica_CompressaoTemSinalNegativo(tc)
            E = 200000;  A = 50;  L = 2000;  P = 10000;
            prob.nodes0   = [0, L; 0, 0];
            prob.elements = [1; 2];
            prob.apoios   = [1 2 4];
            prob.E        = E;
            prob.dens     = 7.85e-6;
            prob.F_total  = [0; 0; -P; 0];   % compressao

            [~, Sigma] = fem_linear_solver(prob, A);
            tc.verifyEqual(Sigma(1), -P/A, 'RelTol', 1e-12, ...
                'Convencao de sinal: compressao deve ser negativa.');
        end

        function testMassaDeBarraUnica(tc)
            E = 200000;  A = 50;  L = 2000;  dens = 7.85e-6;
            prob.nodes0   = [0, L; 0, 0];
            prob.elements = [1; 2];
            prob.apoios   = [1 2 4];
            prob.E        = E;
            prob.dens     = dens;
            prob.F_total  = zeros(4,1);

            w = fem_linear_solver(prob, A);
            tc.verifyEqual(w, A*L*dens, 'RelTol', 1e-12);
        end

        % =================================================================
        % 2. SIMETRIA
        % =================================================================
        function testTrelicaSimetrica_ProduzRespostaSimetrica(tc)
            % Trelica de 2 barras simetrica com carga vertical no apice:
            % as duas barras devem ter tensoes IDENTICAS e o apice nao pode
            % se deslocar horizontalmente.
            caso = problema_trelica_rasa_2barras();
            prob = caso;
            prob.F_total = zeros(6,1);
            prob.F_total(caso.dof_apice_v) = -1000;

            [~, Sigma, ~, d] = fem_linear_solver(prob, caso.areas);

            tc.verifyEqual(Sigma(1), Sigma(2), 'RelTol', 1e-12, ...
                'Barras simetricas deveriam ter a mesma tensao.');
            tc.verifyEqual(d.u(caso.dof_apice_h), 0, 'AbsTol', 1e-9, ...
                'Apice nao deveria se deslocar horizontalmente sob carga simetrica.');
        end

        function testTrelicaSimetrica_CargaVerticalGeraCompressao(tc)
            % Carga para baixo no apice de uma trelica em "A": ambas as
            % barras vao a compressao.
            caso = problema_trelica_rasa_2barras();
            prob = caso;
            prob.F_total = zeros(6,1);
            prob.F_total(caso.dof_apice_v) = -1000;

            [~, Sigma] = fem_linear_solver(prob, caso.areas);
            tc.verifyLessThan(Sigma(1), 0, 'Barra 1 deveria estar comprimida.');
            tc.verifyLessThan(Sigma(2), 0, 'Barra 2 deveria estar comprimida.');
        end

        % =================================================================
        % 3. EQUILIBRIO E LINEARIDADE
        % =================================================================
        function testEquilibrio_KuIgualF(tc)
            % Nos GDLs livres, K*u deve reproduzir o vetor de forcas.
            caso = problema_hadi_10barras();
            [~,~,~,d] = fem_linear_solver(caso, caso.ref_areas);

            residuo = d.K * d.u - caso.F_total;
            residuo(caso.apoios) = 0;

            tc.verifyLessThan(norm(residuo)/norm(caso.F_total), 1e-10, ...
                'K*u nao reproduz F nos GDLs livres.');
        end

        function testLinearidade_DobrarCargaDobraResposta(tc)
            % Propriedade definidora da analise linear.
            caso = problema_hadi_10barras();

            [~, S1, u1] = fem_linear_solver(caso, caso.ref_areas);

            prob = caso;
            prob.F_total = caso.F_total * 2;
            [~, S2, u2] = fem_linear_solver(prob, caso.ref_areas);

            tc.verifyEqual(u2, 2*u1, 'RelTol', 1e-12, ...
                'Dobrar a carga deveria dobrar os deslocamentos.');
            tc.verifyEqual(S2, 2*S1, 'RelTol', 1e-12, ...
                'Dobrar a carga deveria dobrar as tensoes.');
        end

        function testMatrizDeRigidezEhSimetrica(tc)
            caso = problema_hadi_10barras();
            [~,~,~,d] = fem_linear_solver(caso, caso.ref_areas);
            tc.verifyEqual(d.K, d.K', 'AbsTol', 1e-6, ...
                'A matriz de rigidez deve ser simetrica.');
        end

        function testApoiosNaoSeDeslocam(tc)
            caso = problema_hadi_10barras();
            [~,~,~,d] = fem_linear_solver(caso, caso.ref_areas);
            tc.verifyEqual(d.u(caso.apoios), zeros(numel(caso.apoios),1), ...
                'AbsTol', 1e-12, 'GDLs restringidos devem ter deslocamento nulo.');
        end

        % =================================================================
        % 4. VALIDACAO DE ENTRADA
        % =================================================================
        function testRejeitaAreasComTamanhoErrado(tc)
            caso = problema_hadi_10barras();
            tc.verifyError(@() fem_linear_solver(caso, ones(1,3)), ...
                'fem_linear_solver:areasIncompativeis');
        end

        function testRejeitaProblemaIncompleto(tc)
            caso = problema_hadi_10barras();
            caso = rmfield(caso, 'dens');
            tc.verifyError(@() fem_linear_solver(caso, ones(1,10)), ...
                'fem_linear_solver:campoAusente');
        end

        % =================================================================
        % 5. EQUIVALENCIA COM O SOLVER ORIGINAL (arquivado em 05_legado)
        % =================================================================
        function testSolverOriginalProduzResultadoIdentico(tc)
            caso  = problema_hadi_10barras();
            areas = caso.ref_areas;

            [w1, S1, u1] = fem_truss_linear_pre_bloco2(areas);
            [w2, S2, u2] = fem_linear_solver(caso, areas);

            tc.verifyEqual(w1, w2, 'AbsTol', 0);
            tc.verifyEqual(S1, S2, 'AbsTol', 0);
            tc.verifyEqual(u1, u2, 'AbsTol', 0);
        end
    end
end
