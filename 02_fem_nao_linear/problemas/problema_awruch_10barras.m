function caso = problema_awruch_10barras()
% PROBLEMA_AWRUCH_10BARRAS  Trelica de 10 barras com catalogos independentes por barra.
%
% BLOCO 2b do projeto CPIO III — ARQUIVO DE PARAMETROS (sem logica de solver).
%
% =========================================================================
% ORIGEM E STATUS DE VALIDACAO
% =========================================================================
% ATENCAO: este catalogo NAO possui artigo-fonte. Foi montado pelo professor
% da disciplina apenas para EXERCITAR o solver com catalogos distintos por
% barra (em vez de um catalogo unico compartilhado, como no benchmark de
% Hadi & Alvani). Portanto:
%
%   - NAO ha solucao de referencia publicada para comparar.
%   - NAO deve ser usado como evidencia de validacao do solver.
%   - Serve como teste de FLEXIBILIDADE: confirma que o PSO-RID lida com
%     variaveis discretas de catalogos e cardinalidades diferentes entre si.
%
% CORRECAO APLICADA (Bloco 2/3): a versao anterior deste arquivo
% (05_legado/pre_bloco3/catalogo_awruch.m) declarava os campos
%   ref_areas = [19355, 65, 16129, 7742, 65, 65, 5161, 16129, 12903, 65]
%   ref_peso  = 2325.2
% COPIADOS do benchmark de Hadi. Esses valores estavam ERRADOS neste
% contexto: nenhuma daquelas areas pertence aos catalogos definidos abaixo
% (verificado numericamente). Eram codigo morto e enganoso, e foram
% removidos. O peso de referencia fica como NaN, sinalizando "sem referencia".
%
% =========================================================================
% GEOMETRIA E CARREGAMENTO
% =========================================================================
% Identicos aos do benchmark de Hadi & Alvani (2003) — o que muda e apenas o
% conjunto de secoes disponiveis para cada barra. A geometria e reutilizada
% de problema_hadi_10barras para nao duplicar dados.
%
% SAIDA
%   caso : struct pronto para os solvers do Bloco 2a, com .config_vars
%          definindo 10 variaveis discretas de catalogos independentes.

% Reaproveita geometria, material, apoios, cargas e limites do benchmark.
caso = problema_hadi_10barras();

caso.nome = 'Catalogo Awruch - Trelica 10 Barras (catalogos por barra)';

% -------------------------------------------------------------------------
% CATALOGOS POR GRUPO DE BARRAS
% -------------------------------------------------------------------------
% Valores originais fornecidos em polegadas quadradas, convertidos para mm^2.
POL2_PARA_MM2 = 645.16;   % 1 in^2 = 645.16 mm^2 (exato, por definicao da polegada)

cat_A1 = [21.5, 22.5, 23.5, 24.5]             * POL2_PARA_MM2;
cat_A2 = [0.1,  0.15, 0.2,  0.25]             * POL2_PARA_MM2;
cat_A3 = [22.4, 25.4, 27.4, 29.4]             * POL2_PARA_MM2;
cat_A4 = [14.1, 14.2, 14.3, 14.4, 14.5]       * POL2_PARA_MM2;
cat_A6 = [0.5,  1.0,  1.5,  2.0,  2.5]        * POL2_PARA_MM2;
cat_A7 = [11.0, 11.3, 11.7, 12.0, 12.3, 12.5] * POL2_PARA_MM2;

caso.catalogos_por_barra = { ...
    cat_A1, cat_A2, cat_A3, cat_A4, cat_A2, ...
    cat_A6, cat_A7, cat_A3, cat_A4, cat_A2 };

% Remove os campos herdados que NAO se aplicam a este caso
caso = rmfield(caso, 'catalogo');    % nao ha catalogo unico compartilhado

% -------------------------------------------------------------------------
% SEM SOLUCAO DE REFERENCIA (ver nota de origem acima)
% -------------------------------------------------------------------------
caso.ref_areas = [];
caso.ref_peso  = NaN;

% -------------------------------------------------------------------------
% VARIAVEIS DE PROJETO (para o PSO-RID do Bloco 1)
% -------------------------------------------------------------------------
caso.config_vars = struct('tipo', {}, 'opcoes', {});
for i = 1:numel(caso.catalogos_por_barra)
    caso.config_vars(i).tipo   = 'D';
    caso.config_vars(i).opcoes = caso.catalogos_por_barra{i};
end

end
