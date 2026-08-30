function caso = problema_hadi_10barras()
% PROBLEMA_HADI_10BARRAS  Parametros do benchmark de trelica plana de 10 barras.
%
% BLOCO 2b do projeto CPIO III — ARQUIVO DE PARAMETROS (sem logica de solver).
%
% =========================================================================
% REFERENCIA
% =========================================================================
% [HA2003] Hadi, M.N.S.; Alvani, K.S. "Discrete Optimum Design of
%          Geometrically Non-Linear Trusses using Genetic Algorithms".
%          Civil-Comp Press, 2003, Paper 37.  ->  Secao 6.1, Figura 1, Tabela 1
%          -> 00_docs/artigos/Discrete Optimum Design of Geometrically
%             Non-Linear Trusses using Genetic Algorithms.PDF
%
% =========================================================================
% GEOMETRIA (Figura 1 de [HA2003])
% =========================================================================
%
%      |<--- 9144 mm --->|<--- 9144 mm --->|
%
%      5 ------(1)------ 3 ------(2)------ 1     ---
%      |  \              |  \              |      ^
%      |    (7)          |    (9)          |      |
%      |       \         |       \         |    9144 mm
%     (8)        \      (5)        \      (6)     |
%      |           \     |           \     |      v
%      6 ------(3)------ 4 ------(4)------ 2     ---
%                        |                 |
%                        v 445.4 kN        v 445.4 kN
%
%   Nos 5 e 6 sao apoios fixos (engastados nas duas direcoes).
%   Cargas verticais para baixo nos nos 2 e 4.
%
% =========================================================================
% PROPRIEDADES (Secao 6.1 de [HA2003])
% =========================================================================
%   "the modulus of elasticity was 6.89 x 10^4 MPa, and the material density
%    2770 kg/m^3 the allowable displacement was limited to 50.8 mm and the
%    allowable stress to 172.25 MPa in both tension and compression for all
%    members. The problem has 10 design variables."
%
% UNIDADES: mm / N / MPa / kg  (densidade convertida para kg/mm^3)
%
% SAIDA
%   caso : struct com os campos exigidos pelos solvers do Bloco 2a
%          (.nodes0 .elements .apoios .F_total .E .dens) mais os dados de
%          otimizacao (.catalogo .sigma_max .d_max .config_vars .ref_*)

caso.nome = 'Hadi & Alvani (2003) - Trelica 10 Barras (Case 2, discreto)';

% -------------------------------------------------------------------------
% MATERIAL
% -------------------------------------------------------------------------
caso.E    = 6.89e4;      % [MPa]      modulo de elasticidade  [HA2003] Sec. 6.1
caso.dens = 2.770e-6;    % [kg/mm^3]  2770 kg/m^3 = 2.770e-6 kg/mm^3

% -------------------------------------------------------------------------
% GEOMETRIA — coordenadas dos nos [mm]
% -------------------------------------------------------------------------
% O modulo da malha e 9144 mm (= 360 in), conforme a Figura 1 de [HA2003].
caso.modulo = 9144;

%             no    1  2  3  4  5  6
caso.nodes0 = [     2  2  1  1  0  0 ;    % x
                    1  0  1  0  1  0 ] * caso.modulo;   % y

% -------------------------------------------------------------------------
% TOPOLOGIA — conectividade das 10 barras
% -------------------------------------------------------------------------
%   barra     1  2  3  4  5  6  7  8  9 10
caso.elements = [5  3  6  4  3  1  5  6  3  4 ;   % no inicial
                 3  1  4  2  4  2  4  3  2  1];   % no final

% -------------------------------------------------------------------------
% CONDICOES DE CONTORNO — nos 5 e 6 totalmente restringidos
% -------------------------------------------------------------------------
% GDLs: no k ocupa (2k-1) em x e (2k) em y.
%   no 5 -> GDLs  9 e 10
%   no 6 -> GDLs 11 e 12
caso.apoios = [9 10 11 12];

% -------------------------------------------------------------------------
% CARREGAMENTO
% -------------------------------------------------------------------------
% [HA2003] Figura 1 indica 445.4 kN nos nos 2 e 4 (vertical, para baixo).
% Esse valor e o arredondamento de 100 kip:
%   100 kip x 4.4482216152605 N/lbf x 1000 = 444 822.16 N
%
% ATENCAO (discrepancia documentada): o artigo rotula a carga como 445.4 kN,
% mas o valor exato de 100 kip e 444.822 kN — diferenca de 0.13%. Usa-se aqui
% o valor EXATO da conversao. Trocar para 445400 N piora ligeiramente a
% violacao de deslocamento da solucao de referencia. Ver
% 00_docs/notas_e_relatorios/ para a analise dessa discrepancia.
caso.P_no       = 100 * 4.4482216152605 * 1000;   % [N] carga por no carregado
caso.nos_carga  = [2 4];                          % nos 2 e 4

s_dof = 2 * size(caso.nodes0, 2);
caso.F_total = zeros(s_dof, 1);
for k = caso.nos_carga
    caso.F_total(2*k) = -caso.P_no;   % componente y, sentido negativo
end

% -------------------------------------------------------------------------
% RESTRICOES DE PROJETO  ([HA2003] Sec. 6.1)
% -------------------------------------------------------------------------
caso.sigma_max = 172.25;   % [MPa] tensao admissivel (tracao e compressao)
caso.d_max     =  50.80;   % [mm]  deslocamento nodal admissivel

% -------------------------------------------------------------------------
% CATALOGO DE SECOES  ([HA2003] Sec. 6.1, "second catalogue")
% -------------------------------------------------------------------------
% "the second catalogue was: (65, 645, 1290, 3226, 5161, 7742, 9677, 11613,
%  12903, 16129, 19355, 22581, 29032) (mm^2)"
caso.catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, ...
                 11613, 12903, 16129, 19355, 22581, 29032];

% Limites para a variante CONTINUA do problema ([HA2003] Sec. 6.1:
% "the lower and upper bound of cross-sectional areas were given as 65 mm^2
%  and 29032 mm^2")
caso.lb = 65;
caso.ub = 29032;

% -------------------------------------------------------------------------
% SOLUCAO DE REFERENCIA  ([HA2003] Tabela 1, Case 2, coluna "Discrete")
% -------------------------------------------------------------------------
% A1=19355  A2=65     A3=16129  A4=7742   A5=65
% A6=65     A7=5161   A8=16129  A9=12903  A10=65     Peso = 2325.2 kg
caso.ref_areas = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65];
caso.ref_peso  = 2325.2;   % [kg]

% -------------------------------------------------------------------------
% CONFIGURACAO DAS VARIAVEIS DE PROJETO (para o PSO-RID do Bloco 1)
% -------------------------------------------------------------------------
% 10 variaveis discretas, todas escolhendo do mesmo catalogo.
n_barras = size(caso.elements, 2);
for i = 1:n_barras
    caso.config_vars(i).tipo   = 'D';
    caso.config_vars(i).opcoes = caso.catalogo;
end

end
