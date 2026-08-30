classdef TestDecodificadorPSO < matlab.unittest.TestCase
% TESTDECODIFICADORPSO  Testes de rid_mapear_dimensoes e rid_decodificar.
%
% BLOCO TESTADO: 1 (PSO-RID). Nenhum FEM e chamado aqui — os auxiliares sao
% exercitados DIRETAMENTE, permitindo revisao manual peca por peca.
%
% REFERENCIA:
%   [DF2011] Datta, D.; Figueira, J.R. "A real-integer-discrete-coded
%            particle swarm optimization for design problems".
%            Applied Soft Computing 11 (2011) 3625-3633.
%
% Executar:  runtests('06_testes/TestDecodificadorPSO.m')

    methods (TestClassSetup)
        function adicionarCaminhos(~)
            aqui = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(aqui, '..', '01_pso_rid')));
        end
    end

    methods (Test)

        % =================================================================
        % Eq. (5) [DF2011] — dimensionamento da particula
        % =================================================================
        function testEq5_NumeroDeBitsDiscreto(tc)
            % [DF2011] Sec. 4.1: discreta = inteira cujo valor e o INDICE do
            % valor real, com limites 1..N. Logo x_int vai de 0 a N-1 e
            % B = ceil(log2(N)).
            %
            % Conferencia com o artigo: Sec. 5.2 (mola de compressao), o
            % diametro do arame tem 42 valores discretos e usa "6 binary bits".
            casos = { 2,1;  3,2;  4,2;  5,3;  13,4;  16,4;  17,5;  42,6 };
            for i = 1:size(casos,1)
                N = casos{i,1};  esperado = casos{i,2};
                cfg = struct('tipo','D','opcoes', 1:N);
                [mapa, total] = rid_mapear_dimensoes(cfg);
                tc.verifyEqual(mapa(1).n_bits, esperado, ...
                    sprintf('Discreta N=%d deveria usar %d bits.', N, esperado));
                tc.verifyEqual(total, esperado);
            end
        end

        function testEq5_NumeroDeBitsInteiro(tc)
            % Inteira: a Eq. (6) precisa gerar o proprio var.max, logo
            % 2^B - 1 >= max  =>  B = ceil(log2(max+1)).
            %
            % Conferencia com o artigo: Sec. 5.1 (trem de engrenagens),
            % 12 <= z <= 60 usa "6 binary bits" -> ceil(log2(61)) = 6.
            % Sec. 5.2 (mola), N inteiro em [1,70] usa "7 binary bits".
            casos = { 60,6;  70,7;  63,6;  64,7;  1,1 };
            for i = 1:size(casos,1)
                M = casos{i,1};  esperado = casos{i,2};
                cfg = struct('tipo','I','min',0,'max',M);
                mapa = rid_mapear_dimensoes(cfg);
                tc.verifyEqual(mapa(1).n_bits, esperado, ...
                    sprintf('Inteira max=%d deveria usar %d bits.', M, esperado));
            end
        end

        function testEq5_VariavelRealOcupaUmaDimensaoSemBits(tc)
            cfg = struct('tipo','R','min',0,'max',1);
            [mapa, total] = rid_mapear_dimensoes(cfg);
            tc.verifyEqual(mapa(1).n_bits, 0);
            tc.verifyEqual(total, 1);
        end

        function testEq5_SomatorioMisto(tc)
            % L = R + sum B^(I) + sum B^(D)
            %   1 real (1) + inteira max=60 (6 bits) + discreta N=13 (4 bits) = 11
            cfg(1) = struct('tipo','R','min',0,'max',1,'opcoes',[]);
            cfg(2) = struct('tipo','I','min',12,'max',60,'opcoes',[]);
            cfg(3) = struct('tipo','D','min',[],'max',[],'opcoes',1:13);
            [mapa, total] = rid_mapear_dimensoes(cfg);
            tc.verifyEqual([mapa.n_bits], [0 6 4]);
            tc.verifyEqual(total, 11);
            % Blocos contiguos e sem sobreposicao
            tc.verifyEqual([mapa.inicio], [1 2 8]);
            tc.verifyEqual([mapa.fim],    [1 7 11]);
        end

        % =================================================================
        % Eq. (6) [DF2011] — decodificacao binario -> inteiro
        % =================================================================
        function testEq6_SomaPonderadaExata(tc)
            % x = sum_{i=1}^{B} 2^(B-i) * b_i, bit mais significativo a esquerda.
            % Variavel inteira com folga suficiente para nao violar limites.
            cfg  = struct('tipo','I','min',0,'max',255);
            mapa = rid_mapear_dimensoes(cfg);          % 8 bits
            tc.verifyEqual(mapa(1).n_bits, 8);

            casos = { [0 0 0 0 0 0 0 0],   0 ;
                      [0 0 0 0 0 0 0 1],   1 ;
                      [0 0 0 0 0 0 1 0],   2 ;
                      [1 0 0 0 0 0 0 0], 128 ;
                      [1 0 1 0 1 0 1 0], 170 ;
                      [1 1 1 1 1 1 1 1], 255 };

            for i = 1:size(casos,1)
                bits     = casos{i,1};
                esperado = casos{i,2};
                [v, viol] = rid_decodificar(bits, mapa, cfg);
                tc.verifyEqual(v, esperado, ...
                    sprintf('Eq. (6) falhou para bits [%s].', num2str(bits)));
                tc.verifyEqual(viol, 0, 'Nao deveria haver violacao dentro de [0,255].');
            end
        end

        % =================================================================
        % [D1] Decodificacao discreta — modo 'datta' (fiel ao artigo)
        % =================================================================
        function testD1_Datta_IndiceEhValorMaisUm(tc)
            % [DF2011] Sec. 4.1: o inteiro decodificado e o INDICE do valor
            % discreto, com limite inferior 1. Logo indice = x_int + 1.
            catalogo = [10 20 30 40 50 60 70 80];      % N=8 -> 3 bits, sem overflow
            cfg  = struct('tipo','D','opcoes',catalogo);
            mapa = rid_mapear_dimensoes(cfg);
            tc.verifyEqual(mapa(1).n_bits, 3);

            for x_int = 0:7
                bits = tc.paraBits(x_int, 3);
                [v, viol] = rid_decodificar(bits, mapa, cfg, 'datta');
                tc.verifyEqual(v, catalogo(x_int + 1), ...
                    sprintf('x_int=%d deveria mapear para opcoes(%d).', x_int, x_int+1));
                tc.verifyEqual(viol, 0);
            end
        end

        function testD1_Datta_OverflowGeraViolacaoEstrutural(tc)
            % Catalogo do benchmark Hadi: N=13 -> 4 bits -> 16 estados.
            % Indices 14, 15 e 16 (x_int = 13, 14, 15) NAO existem e devem
            % gerar violacao estrutural proporcional ao excedente.
            catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, ...
                        11613, 12903, 16129, 19355, 22581, 29032];
            cfg  = struct('tipo','D','opcoes',catalogo);
            mapa = rid_mapear_dimensoes(cfg);
            tc.verifyEqual(mapa(1).n_bits, 4);

            % Dentro do catalogo: sem violacao
            for x_int = 0:12
                [v, viol] = rid_decodificar(tc.paraBits(x_int,4), mapa, cfg, 'datta');
                tc.verifyEqual(viol, 0);
                tc.verifyEqual(v, catalogo(x_int+1));
            end

            % Fora do catalogo: violacao = indice - N
            esperadas = [1 2 3];   % x_int 13,14,15 -> idx 14,15,16 -> viol 1,2,3
            for j = 1:3
                x_int = 12 + j;
                [~, viol] = rid_decodificar(tc.paraBits(x_int,4), mapa, cfg, 'datta');
                tc.verifyEqual(viol, esperadas(j), ...
                    sprintf('x_int=%d deveria gerar violacao estrutural %d.', ...
                            x_int, esperadas(j)));
            end
        end

        function testD1_Datta_NuncaVazaValorForaDoCatalogo(tc)
            % Mesmo em overflow, o valor devolvido deve pertencer ao catalogo
            % (e um preenchimento) — o FEM jamais deve receber area inexistente.
            catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, ...
                        11613, 12903, 16129, 19355, 22581, 29032];
            cfg  = struct('tipo','D','opcoes',catalogo);
            mapa = rid_mapear_dimensoes(cfg);

            for x_int = 0:15
                v = rid_decodificar(tc.paraBits(x_int,4), mapa, cfg, 'datta');
                tc.verifyTrue(ismember(v, catalogo), ...
                    sprintf('x_int=%d vazou valor fora do catalogo.', x_int));
            end
        end

        % =================================================================
        % [D1] Decodificacao discreta — modo 'proporcional' (legado)
        % =================================================================
        function testD1_Proporcional_NuncaGeraOverflow(tc)
            catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, ...
                        11613, 12903, 16129, 19355, 22581, 29032];
            cfg  = struct('tipo','D','opcoes',catalogo);
            mapa = rid_mapear_dimensoes(cfg);

            for x_int = 0:15
                [v, viol] = rid_decodificar(tc.paraBits(x_int,4), mapa, cfg, 'proporcional');
                tc.verifyEqual(viol, 0, ...
                    'Modo proporcional nao deveria gerar violacao estrutural.');
                tc.verifyTrue(ismember(v, catalogo));
            end
        end

        function testD1_CaracterizacaoDoTradeOff(tc)
            % TESTE DE CARACTERIZACAO (documenta o trade-off, nao reprova).
            %
            % Quantifica exatamente o custo de cada estrategia sobre o
            % catalogo do benchmark Hadi (N=13, B=4, 16 estados).
            catalogo = 1:13;
            cfg  = struct('tipo','D','opcoes',catalogo);
            mapa = rid_mapear_dimensoes(cfg);
            N = 13;  B = 4;  n_estados = 2^B;

            % --- modo datta: conta codigos invalidos e mede uniformidade ---
            cont_datta = zeros(1, N);
            n_invalidos = 0;
            for x_int = 0:n_estados-1
                [v, viol] = rid_decodificar(tc.paraBits(x_int,B), mapa, cfg, 'datta');
                if viol > 0
                    n_invalidos = n_invalidos + 1;
                else
                    cont_datta(v) = cont_datta(v) + 1;
                end
            end

            % --- modo proporcional: conta ocorrencias por opcao ---
            cont_prop = zeros(1, N);
            for x_int = 0:n_estados-1
                v = rid_decodificar(tc.paraBits(x_int,B), mapa, cfg, 'proporcional');
                cont_prop(v) = cont_prop(v) + 1;
            end

            % datta: entre os codigos VALIDOS a distribuicao e perfeitamente
            % uniforme (cada opcao alcancada por exatamente 1 codigo)
            tc.verifyEqual(cont_datta, ones(1,N), ...
                'No modo datta cada opcao valida deveria ter exatamente 1 codigo.');
            tc.verifyEqual(n_invalidos, n_estados - N);

            % proporcional: sem desperdicio, mas com vies
            tc.verifyEqual(sum(cont_prop), n_estados, ...
                'Modo proporcional deveria usar todos os codigos.');
            vies_prop = max(cont_prop) / min(cont_prop);

            p_valida_10 = (N / n_estados)^10;

            fprintf(['\n[CARACTERIZACAO D1] catalogo N=%d, B=%d bits, %d estados\n' ...
                     '  datta        : %d/%d codigos invalidos (%.2f%%), vies entre validos = 1.00x\n' ...
                     '  proporcional : 0 codigos invalidos, vies max/min = %.2fx\n' ...
                     '  P(solucao de 10 variaveis 100%% valida | datta) = %.2f%%\n'], ...
                    N, B, n_estados, n_invalidos, n_estados, ...
                    100*n_invalidos/n_estados, vies_prop, 100*p_valida_10);

            % Asserts de caracterizacao (travam os numeros documentados)
            tc.verifyEqual(n_invalidos, 3);
            tc.verifyEqual(vies_prop, 2);
            tc.verifyLessThan(p_valida_10, 0.13);
        end

        % =================================================================
        % [D2] Variavel inteira — limites viram violacao, sem clamp
        % =================================================================
        function testD2_RespeitaLimiteInferior(tc)
            % REGRESSAO: a versao anterior fazia min(x, max) e IGNORAVA
            % completamente var.min, deixando passar valores abaixo do limite.
            % [DF2011] Eq. (12): trem de engrenagens com 12 <= z <= 60.
            cfg  = struct('tipo','I','min',12,'max',60);
            mapa = rid_mapear_dimensoes(cfg);          % ceil(log2(61)) = 6 bits
            tc.verifyEqual(mapa(1).n_bits, 6);

            % x_int = 5 esta abaixo de min=12 -> violacao = 12 - 5 = 7
            [v, viol] = rid_decodificar(tc.paraBits(5,6), mapa, cfg);
            tc.verifyEqual(v, 5, 'Valor deve ser reportado sem clamp (Eq. 6).');
            tc.verifyEqual(viol, 7, 'Violacao do limite INFERIOR nao detectada (bug D2).');
        end

        function testD2_RespeitaLimiteSuperior(tc)
            cfg  = struct('tipo','I','min',12,'max',60);
            mapa = rid_mapear_dimensoes(cfg);

            % x_int = 63 excede max=60 -> violacao = 3
            [v, viol] = rid_decodificar(tc.paraBits(63,6), mapa, cfg);
            tc.verifyEqual(v, 63, 'Valor deve ser reportado sem clamp (Eq. 6).');
            tc.verifyEqual(viol, 3);
        end

        function testD2_DentroDosLimitesNaoViola(tc)
            cfg  = struct('tipo','I','min',12,'max',60);
            mapa = rid_mapear_dimensoes(cfg);
            for x_int = 12:60
                [v, viol] = rid_decodificar(tc.paraBits(x_int,6), mapa, cfg);
                tc.verifyEqual(v, x_int);
                tc.verifyEqual(viol, 0);
            end
        end

        % =================================================================
        % Variavel real e casos mistos
        % =================================================================
        function testReal_ValorPassaDireto(tc)
            cfg  = struct('tipo','R','min',-3.5,'max',7.25);
            mapa = rid_mapear_dimensoes(cfg);
            [v, viol] = rid_decodificar(2.7182818, mapa, cfg);
            tc.verifyEqual(v, 2.7182818, 'AbsTol', 1e-12);
            tc.verifyEqual(viol, 0);
        end

        function testMisto_ViolacoesSomam(tc)
            % Duas variaveis fora do dominio: as violacoes devem SOMAR.
            cfg(1) = struct('tipo','I','min',12,'max',60,'opcoes',[]);
            cfg(2) = struct('tipo','D','min',[],'max',[],'opcoes',1:13);
            mapa = rid_mapear_dimensoes(cfg);          % 6 bits + 4 bits

            vetor = [tc.paraBits(63,6), tc.paraBits(15,4)];
            % inteira: 63 > 60      -> violacao 3
            % discreta: idx 16 > 13 -> violacao 3
            [~, viol] = rid_decodificar(vetor, mapa, cfg, 'datta');
            tc.verifyEqual(viol, 6, 'Violacoes de variaveis distintas devem somar.');
        end

        function testRejeitaModoInvalido(tc)
            cfg  = struct('tipo','D','opcoes',1:4);
            mapa = rid_mapear_dimensoes(cfg);
            tc.verifyError(@() rid_decodificar([0 0], mapa, cfg, 'inexistente'), ...
                'rid_decodificar:modoInvalido');
        end
    end

    % =====================================================================
    methods (Access = private)
        function bits = paraBits(~, valor, n_bits)
            % Converte inteiro para vetor de bits, mais significativo a
            % esquerda (mesma convencao da Eq. 6 de [DF2011]).
            bits = zeros(1, n_bits);
            for i = n_bits:-1:1
                bits(i) = mod(valor, 2);
                valor   = floor(valor / 2);
            end
        end
    end
end
