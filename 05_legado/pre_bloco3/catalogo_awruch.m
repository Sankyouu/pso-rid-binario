function caso = catalogo_awruch()
% CATALOGO_AWRUCH Retorna catálogos por grupos de barras de Awruch
% Treliça de 10 barras com seções discretas por elemento

caso.nome = 'Catálogo Awruch - Treliça 10 Barras';

in2_to_mm2 = 645.16; % Fator de conversão

cat_A1 = [21.5, 22.5, 23.5, 24.5]             * in2_to_mm2;
cat_A2 = [0.1,  0.15, 0.2,  0.25]             * in2_to_mm2;
cat_A3 = [22.4, 25.4, 27.4, 29.4]             * in2_to_mm2;
cat_A4 = [14.1, 14.2, 14.3, 14.4, 14.5]       * in2_to_mm2;
cat_A6 = [0.5,  1.0,  1.5,  2.0,  2.5]        * in2_to_mm2;
cat_A7 = [11.0, 11.3, 11.7, 12.0, 12.3, 12.5] * in2_to_mm2;

caso.sigma_max = 172.25; % [MPa]
caso.d_max     =  50.80; % [mm]
caso.ref_peso  =    NaN;

% Configuração das variáveis discretas por barra
caso.config_vars(1).tipo   = 'D'; caso.config_vars(1).opcoes = cat_A1;
caso.config_vars(2).tipo   = 'D'; caso.config_vars(2).opcoes = cat_A2;
caso.config_vars(3).tipo   = 'D'; caso.config_vars(3).opcoes = cat_A3;
caso.config_vars(4).tipo   = 'D'; caso.config_vars(4).opcoes = cat_A4;
caso.config_vars(5).tipo   = 'D'; caso.config_vars(5).opcoes = cat_A2;
caso.config_vars(6).tipo   = 'D'; caso.config_vars(6).opcoes = cat_A6;
caso.config_vars(7).tipo   = 'D'; caso.config_vars(7).opcoes = cat_A7;
caso.config_vars(8).tipo   = 'D'; caso.config_vars(8).opcoes = cat_A3;
caso.config_vars(9).tipo   = 'D'; caso.config_vars(9).opcoes = cat_A4;
caso.config_vars(10).tipo  = 'D'; caso.config_vars(10).opcoes = cat_A2;

% Solução de referência exata da Tabela 1 de Hadi & Alvani (2003), Case 2 Discreto:
% Áreas [mm²]: A1=19355, A2=65, A3=16129, A4=7742, A5=65, A6=65, A7=5161, A8=16129, A9=12903, A10=65
caso.ref_areas = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65];
caso.ref_peso  = 2325.2; % [kg] (conforme reportado na Tabela 1 do artigo)

end
