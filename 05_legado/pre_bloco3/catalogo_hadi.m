function caso = catalogo_hadi()
% CATALOGO_HADI Retorna parâmetros e catálogo do benchmark de Hadi & Alvani (2003)
% Artigo: "Discrete Optimum Design of Geometrically Non-Linear Trusses using Genetic Algorithms"
% Treliça de 10 barras, Case 2 (Variáveis Discretas)

caso.nome = 'Hadi & Alvani (2003) - Treliça 10 Barras (Case 2)';

% Catálogo de 13 seções comerciais disponíveis [mm²]
caso.catalogo = [65, 645, 1290, 3226, 5161, 7742, 9677, 11613, 12903, 16129, 19355, 22581, 29032];

% Limites das restrições (Seção 6.1 do artigo)
caso.sigma_max = 172.25; % [MPa] — tensão admissível tração/compressão
caso.d_max     =  50.80; % [mm]  — deslocamento nodal máximo admissível

% Limites contínuos (para testes contínuos)
caso.lb = 65;
caso.ub = 29032;

% Solução de referência exata da Tabela 1 de Hadi & Alvani (2003), Case 2 Discreto:
% Áreas [mm²]: A1=19355, A2=65, A3=16129, A4=7742, A5=65, A6=65, A7=5161, A8=16129, A9=12903, A10=65
caso.ref_areas = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65];
caso.ref_peso  = 2325.2; % [kg] (conforme reportado na Tabela 1 do artigo)

% Configuração das 10 variáveis discretas para o PSO-RID
n_barras = 10;
for i = 1:n_barras
    caso.config_vars(i).tipo   = 'D';
    caso.config_vars(i).opcoes = caso.catalogo;
end

end
