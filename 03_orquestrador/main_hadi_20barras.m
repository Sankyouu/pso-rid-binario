function resultado = main_hadi_20barras(seed, n_runs, modo)
% MAIN_HADI_20BARRAS  Otimizacao da trelica plana de 20 barras de [HA2003] Sec. 2.2.
%
% -------------------------------------------------------------------------
% O QUE ESTE ARQUIVO FAZ
% -------------------------------------------------------------------------
%
% Segundo exemplo de projeto de [HA2003]. Amarra os mesmos blocos do caso de
% 10 barras:
%
%   Bloco 1 (otimizador) : pso_rid.m                <- 01_pso_rid/
%   Bloco 2  (solver)    : fem_nao_linear_solver.m  <- 02_fem_nao_linear/
%                          fem_linear_solver.m
%
% O QUE ESTE CASO TRAZ DE NOVO em relacao aos outros orquestradores:
%
%   1. AGRUPAMENTO DE BARRAS. As 20 barras sao divididas em 10 grupos que
%      compartilham a mesma area. Logo ha 10 variaveis de projeto para 20
%      barras — o PSO otimiza no espaco dos GRUPOS e a funcao objetivo
%      expande para as barras antes de chamar o FEM (ver expandir_areas).
%      Nenhum outro caso do projeto tinha isso: em Hadi-10 e Awruch a
%      correspondencia barra <-> variavel era 1:1.
%
%   2. CATALOGO DE 65 SECOES (contra 13 em Hadi-10). Isso muda o
%      dimensionamento binario: ceil(log2(65)) = 7 bits por variavel, contra
%      4 antes. A particula passa a ter 10 x 7 = 70 dimensoes.
%
%   3. TRELICA BIAPOIADA. Hadi-10 e Awruch sao balanços engastados; esta e
%      simplesmente apoiada (apoio fixo num extremo, movel no outro).
%
% FONTE UNICA DO CASO: a definicao do benchmark vive na funcao local
% caso_hadi_20barras, exposta por main_hadi_20barras('caso') — mesmo padrao
% de main_hadi_nao_linear e main_awruch_discreto.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%
%   [HA2003]  Hadi & Alvani (2003), Civil-Comp Press, Paper 37.
%             Sec. 2.2 "Design of 20-bar Plane Truss", Figura 2, Tabela 2.
%             NOTA DE NUMERACAO: a secao e mesmo "2.2" no PDF, embora o
%             exemplo anterior (10 barras) seja "6.1" e o seguinte (51
%             barras) "6.3". A numeracao do artigo publicado esta furada
%             nesse ponto; o titulo literal e "2.2 Design of 20-bar Plane
%             Truss" e fica entre 6.1 e 6.3.
%   [DF2011]  Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_hadi_20barras                      % semente 42, 5 runs, NAO LINEAR
%   main_hadi_20barras(123)                 % outra semente
%   main_hadi_20barras(42, 10)              % 10 execucoes
%   main_hadi_20barras(42, 5, 'linear')     % analise linear (Case 2)
%   caso = main_hadi_20barras('caso');      % so o struct do problema
%   r = main_hadi_20barras;                 % devolve struct com resultados
%
% SAIDA
%   resultado : struct com .melhor_areas (10 grupos) .melhor_areas_barras (20)
%               .melhor_peso .violacao .Sigma .u .historicos .pesos_por_run
%               .caso .modo

% -------------------------------------------------------------------------
% ACESSOR DO CASO (ver cabecalho)
% -------------------------------------------------------------------------

if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    resultado = caso_hadi_20barras();
    return;
end

% Acesso as funcoes locais para teste, no mesmo padrao de pso_rid('auxiliares').
% Existe para que 06_testes/ possa verificar a reproducao dos valores
% publicados sem reimplementar a funcao objetivo aqui.
if nargin == 1 && ischar(seed) && strcmp(seed, 'auxiliares')
    resultado = struct( ...
        'expandir_areas',  @expandir_areas,  ...
        'avaliar_projeto', @avaliar_projeto);
    return;
end

if nargin < 1 || isempty(seed),   seed   = 42;            end
if nargin < 2 || isempty(n_runs), n_runs = 5;             end
if nargin < 3 || isempty(modo),   modo   = 'nao_linear';  end

assert(any(strcmp(modo, {'linear', 'nao_linear'})), ...
    'main_hadi_20barras:modoInvalido', ...
    'modo deve ser ''linear'' ou ''nao_linear''. Recebido: ''%s''.', modo);

garantir_caminhos();
rng(seed);

% -------------------------------------------------------------------------
% 1. PROBLEMA
% -------------------------------------------------------------------------

caso = caso_hadi_20barras();

% Referencia da Tabela 2 correspondente ao modo de analise escolhido:
%   Case 1 = solucao discreta NAO LINEAR    -> 6222 kg
%   Case 2 = solucao discreta LINEAR        -> 6253 kg
if strcmp(modo, 'nao_linear')
    ref_areas = caso.ref_areas_nao_linear;
    ref_peso  = caso.ref_peso_nao_linear;
else
    ref_areas = caso.ref_areas_linear;
    ref_peso  = caso.ref_peso_linear;
end

% -------------------------------------------------------------------------
% 2. FUNCAO OBJETIVO
% -------------------------------------------------------------------------

funcao_objetivo = @(areas_grupos) avaliar_projeto(areas_grupos, caso, modo);

% -------------------------------------------------------------------------
% 3. CONFIGURACAO DO PSO-RID (Bloco 1)
% -------------------------------------------------------------------------
% Mesmos valores de [DF2011] Sec. 5 usados nos demais orquestradores. O
% espaco de busca aqui e MAIOR que o de Hadi-10 (70 dimensoes binarias
% contra 40, 65^10 combinacoes contra 13^10), entao o orcamento de
% iteracoes e o mesmo mas a estagnacao e mais folgada.
% -------------------------------------------------------------------------

pso_params = struct();
pso_params.n_particulas           = 100;
pso_params.max_iter               = 1000;
pso_params.w                      = 1.0;      % [DF2011] Sec. 5
pso_params.c1                     = 1.0;      % [DF2011] Sec. 5
pso_params.c2                     = 2.0;      % [DF2011] Sec. 5
pso_params.pm                     = 0.15;     % [DF2011] Sec. 4.3/5
pso_params.auto_adaptativo        = true;     % [DF2011] Sec. 5
pso_params.decodificacao_discreta = 'proporcional';  % ver [D1] em pso_rid.m
pso_params.tol_estagnacao         = 250;
pso_params.verbose                = true;
pso_params.print_interval         = 100;

% -------------------------------------------------------------------------
% 4. EXECUCAO MULTI-START
% -------------------------------------------------------------------------

imprimir_cabecalho(caso, seed, n_runs, pso_params, modo, ref_peso);

historicos    = cell(n_runs, 1);
pesos_por_run = nan(n_runs, 1);
melhor_peso   = inf;
melhor_areas  = [];

for r = 1:n_runs
    fprintf('--- Executando Run %d/%d ---\n', r, n_runs);

    [areas_r, peso_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}    = hist_r;
    pesos_por_run(r) = peso_r;

    % Seleciona o melhor apenas entre solucoes VIAVEIS ([DEB2000] criterio 1)
    if det_r.gbest_viol <= 0 && peso_r < melhor_peso
        melhor_peso  = peso_r;
        melhor_areas = areas_r;
    end

    fprintf('    -> Run %d: peso = %.2f kg | violacao = %.3e | iters = %d | avaliacoes FEM = %d\n\n', ...
            r, peso_r, det_r.gbest_viol, det_r.iter_executadas, det_r.n_avaliacoes);
end

if isempty(melhor_areas)
    error('main_hadi_20barras:semSolucaoViavel', ...
        ['Nenhuma das %d execucoes encontrou solucao VIAVEL. ' ...
         'Aumente max_iter/n_particulas ou revise as restricoes.'], n_runs);
end

% -------------------------------------------------------------------------
% 5. AVALIACAO DETALHADA DA MELHOR SOLUCAO
% -------------------------------------------------------------------------

melhor_areas_barras = expandir_areas(melhor_areas, caso);

if strcmp(modo, 'nao_linear')
    [peso_final, Sigma_final, u_final] = fem_nao_linear_solver(caso, melhor_areas_barras);
else
    [peso_final, Sigma_final, u_final] = fem_linear_solver(caso, melhor_areas_barras);
end
[~, viol_final] = avaliar_projeto(melhor_areas, caso, modo);

% -------------------------------------------------------------------------
% 6. RELATORIO E GRAFICOS
% -------------------------------------------------------------------------
% ATENCAO na leitura do relatorio: a coluna rotulada "Barra" traz aqui os 10
% GRUPOS, nao as 20 barras (ver o mapa de grupos impresso logo abaixo).
% -------------------------------------------------------------------------

relatorio_comparativo(sprintf('Trelica 20 Barras Hadi (%s)', upper(strrep(modo,'_',' '))), ...
                      melhor_areas, peso_final, ...
                      ref_areas, ref_peso, ...
                      Sigma_final, caso.sigma_max, ...
                      u_final, caso.d_max, viol_final);

imprimir_mapa_grupos(melhor_areas, caso);
imprimir_estatisticas_runs(pesos_por_run, ref_peso);
verificar_areas_pertencem_ao_catalogo(melhor_areas, caso);

fig = plot_convergencia(historicos, ref_peso, ...
                        sprintf('PSO-RID — Trelica 20 Barras (%s)', ...
                                upper(strrep(modo,'_',' '))));
salvar_figura(fig, sprintf('convergencia_hadi_20barras_%s', modo));

% -------------------------------------------------------------------------
% 7. SAIDA ESTRUTURADA
% -------------------------------------------------------------------------

resultado.caso                = caso;
resultado.modo                = modo;
resultado.melhor_areas        = melhor_areas;         % 1 x 10 (grupos)
resultado.melhor_areas_barras = melhor_areas_barras;  % 1 x 20 (barras)
resultado.melhor_peso         = peso_final;
resultado.violacao            = viol_final;
resultado.Sigma               = Sigma_final;
resultado.u                   = u_final;
resultado.historicos          = historicos;
resultado.pesos_por_run       = pesos_por_run;
resultado.seed                = seed;

end


% #########################################################################
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################


% -------------------------------------------------------------------------
% >>> LOCAL caso_hadi_20barras — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_hadi_20barras('caso')
% -------------------------------------------------------------------------

function caso = caso_hadi_20barras()
% CASO_HADI_20BARRAS  Parametros do benchmark de trelica plana de 20 barras.
% Sem logica de solver.
%
% -------------------------------------------------------------------------
% REFERENCIA
% -------------------------------------------------------------------------
%
% [HA2003] Hadi, M.N.S.; Alvani, K.S. "Discrete Optimum Design of
%          Geometrically Non-Linear Trusses using Genetic Algorithms".
%          Civil-Comp Press, 2003, Paper 37.  ->  Sec. 2.2, Figura 2, Tabela 2
%          -> 00_docs/artigos/Discrete Optimum Design of Geometrically
%             Non-Linear Trusses using Genetic Algorithms.PDF  (pagina 11)
%
% -------------------------------------------------------------------------
% GEOMETRIA (Figura 2 de [HA2003])
% -------------------------------------------------------------------------
%
%   Vao total 5 @ 5000 mm = 25000 mm.  Altura 5000 mm.
%   Banzo inferior: nos 1-2-4-6-8-10.  Banzo superior: nos 3-5-7-9.
%   Os nos superiores ficam exatamente acima dos nos 2, 4, 6 e 8.
%
%          3 ---(5)--- 5 ---(10)--- 7 ---(15)--- 9
%         /|\         /|\          /|\          /|\
%        / | \ (6)   / | \ (11)   / | \ (16)   / | \
%     (2)  |  \     /  |  \      /  |  \      /  |  \ (20)
%      /  (3)  \(7)/  (9)  \(12)/ (14)  \(17)/ (19)  \
%     /    |    \ /    |    \  /    |    \  /    |    \
%    1---(1)----2--(4)-4--(8)-6--(13)-8--(18)----10
%    ^               |         |                  o
%   fixo             v         v                movel
%                 850 kN    850 kN
%
%   Diagonais cruzadas em cada um dos tres paineis internos:
%     painel 2-3-5-4 : (6) = 3->4  e  (7) = 2->5
%     painel 4-5-7-6 : (11) = 5->6 e  (12) = 4->7
%     painel 6-7-9-8 : (16) = 7->8 e  (17) = 6->9
%
%   As demais barras: banzo inferior (1)(4)(8)(13)(18), banzo superior
%   (5)(10)(15), montantes (3)(9)(14)(19), diagonais de extremidade (2)(20).
%
% COMO A TOPOLOGIA FOI CONFIRMADA. A Figura 2 e um desenho, nao uma tabela
% de coordenadas — a conectividade foi lida da figura (pagina 11 do PDF) e
% depois VERIFICADA numericamente: com as coordenadas e o agrupamento abaixo,
% o peso das duas colunas da Tabela 2 e reproduzido exatamente,
%
%     Case 1 (nao linear) : 6222.0 kg calculado  vs  6222 kg publicado
%     Case 2 (linear)     : 6253.0 kg calculado  vs  6253 kg publicado
%
% o que so acontece se comprimentos e grupos estiverem certos (o peso e
% dens * sum(A_i * L_i), logo qualquer erro de coordenada ou de agrupamento
% apareceria na terceira casa). Reproduza com:
%
%   c = main_hadi_20barras('caso');
%   L = arrayfun(@(e) norm(c.nodes0(:,c.elements(2,e)) - ...
%                          c.nodes0(:,c.elements(1,e))), 1:20);
%   a = c.ref_areas_nao_linear(c.grupo_por_barra);
%   c.dens * sum(a .* L)        % -> 6222.145 kg
%
% -------------------------------------------------------------------------
% PROPRIEDADES (Sec. 2.2 de [HA2003])
% -------------------------------------------------------------------------
%   "The modulus of elasticity is specified as 210x10^3 MPa, and material
%    density as 7850 kg/m^3. The allowable displacement is limited to 50 mm
%    and the yield stress according to Ref. [19] for the circular tubeline
%    steel used is 350 MPa."
%
% Ao contrario do caso de 10 barras, aqui NAO ha conversao imperial
% escondida: 210 GPa, 7850 kg/m^3, 50 mm e 350 MPa sao valores metricos
% redondos de aco estrutural. Nao ha discrepancia a documentar.
%
% UNIDADES: mm / N / MPa / kg  (densidade em kg/mm^3)
%
% SAIDA
%   caso : struct com os campos exigidos pelos solvers do Bloco 2
%          (.nodes0 .elements .apoios .F_total .E .dens) mais os dados de
%          otimizacao (.catalogo .sigma_max .d_max .config_vars .grupos
%          .grupo_por_barra .ref_*)

caso.nome = 'Hadi & Alvani (2003) - Trelica 20 Barras (Sec. 2.2)';

% -------------------------------------------------------------------------
% MATERIAL
% -------------------------------------------------------------------------

caso.E    = 210e3;       % [MPa]      210 x 10^3 MPa
caso.dens = 7.850e-6;    % [kg/mm^3]  7850 kg/m^3

% -------------------------------------------------------------------------
% RESTRICOES
% -------------------------------------------------------------------------

caso.sigma_max = 350;    % [MPa] tensao de escoamento (tracao e compressao)
caso.d_max     = 50;     % [mm]  deslocamento admissivel

% -------------------------------------------------------------------------
% GEOMETRIA — coordenadas dos nos [mm]
% -------------------------------------------------------------------------

caso.modulo = 5000;      % [mm] modulo da malha (vao 5@5000, altura 5000)

%             no    1  2  3  4  5  6  7  8  9 10
caso.nodes0 = [     0  1  1  2  2  3  3  4  4  5 ;    % x
                    0  0  1  0  1  0  1  0  1  0 ] * caso.modulo;   % y

% -------------------------------------------------------------------------
% TOPOLOGIA — conectividade das 20 barras
% -------------------------------------------------------------------------

%   barra          1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
caso.elements = [  1  1  2  2  3  3  2  4  4  5  5  4  6  6  7  7  6  8  8  9 ;   % no inicial
                   2  3  3  4  5  4  5  6  5  7  6  7  8  7  9  8  9 10  9 10];   % no final

% -------------------------------------------------------------------------
% CONDICOES DE CONTORNO (Figura 2: apoio fixo no no 1, movel no no 10)
% -------------------------------------------------------------------------

% GDLs: no k ocupa (2k-1) em x e (2k) em y.
%   no  1 -> GDLs  1 e  2   apoio FIXO   (x e y restringidos)
%   no 10 -> GDLs 19 e 20   apoio MOVEL  (so y restringido; livre em x)
caso.apoios = [1 2 20];

% -------------------------------------------------------------------------
% CARREGAMENTO
% -------------------------------------------------------------------------
%
% [HA2003] Sec. 2.2: "The structure is subjected to two load conditions of
% 850 kN each at nodes 4 and 6."
%
% LEITURA ADOTADA: a Figura 2 desenha as DUAS setas simultaneamente, logo
% trata-se de um unico caso de carregamento com duas cargas concentradas de
% 850 kN, e nao de dois casos de carga alternativos. O texto diz "two load
% conditions" mas a figura e inequivoca.

caso.P_no      = 850e3;    % [N] 850 kN por no carregado
caso.nos_carga = [4 6];

s_dof = 2 * size(caso.nodes0, 2);
caso.F_total = zeros(s_dof, 1);
for k = 1:numel(caso.nos_carga)
    caso.F_total(2 * caso.nos_carga(k)) = -caso.P_no;   % vertical, para baixo
end

% -------------------------------------------------------------------------
% AGRUPAMENTO DE BARRAS (Tabela 2 de [HA2003])
% -------------------------------------------------------------------------
%
% As 20 barras sao divididas em 10 grupos de mesma area. O agrupamento
% publicado tem significado estrutural claro, o que serviu de conferencia
% cruzada da topologia lida da figura:
%
%   A1 -> as 5 barras do banzo inferior
%   A2 -> os 4 montantes verticais
%   A3 -> as 3 barras do banzo superior
%   A4 -> as 2 diagonais de extremidade
%   A5..A10 -> as 6 diagonais dos paineis internos, uma por grupo
%
caso.grupos = { ...
    [ 1  4  8 13 18], ...   % A1  banzo inferior
    [ 3  9 14 19],    ...   % A2  montantes
    [ 5 10 15],       ...   % A3  banzo superior
    [ 2 20],          ...   % A4  diagonais de extremidade
    6,  ...                 % A5
    7,  ...                 % A6
    11, ...                 % A7
    12, ...                 % A8
    16, ...                 % A9
    17};                    % A10

caso.nomes_grupos = { ...
    'banzo inferior', 'montantes', 'banzo superior', 'diag. extremidade', ...
    'diag. 3-4', 'diag. 2-5', 'diag. 5-6', 'diag. 4-7', 'diag. 7-8', 'diag. 6-9'};

% Mapa inverso 1 x 20: em que grupo cada barra esta. E o que expandir_areas
% usa — indexacao direta, sem laco sobre grupos a cada avaliacao do FEM.
n_barras = size(caso.elements, 2);
caso.grupo_por_barra = zeros(1, n_barras);
for g = 1:numel(caso.grupos)
    caso.grupo_por_barra(caso.grupos{g}) = g;
end
assert(all(caso.grupo_por_barra > 0), ...
    'caso_hadi_20barras:barraSemGrupo', ...
    'Toda barra precisa pertencer a exatamente um grupo.');

% -------------------------------------------------------------------------
% CATALOGO DE SECOES (Sec. 2.2 de [HA2003], 65 secoes, de Ref. [19])
% -------------------------------------------------------------------------
%
% Transcrito LITERALMENTE do artigo, inclusive na ordem. Note que a lista
% publicada nao e monotonica: "2280, 2230" aparece fora de ordem (posicoes
% 43 e 44). Mantem-se como impresso por fidelidade — o decodificador nao
% exige catalogo ordenado, e apenas a vizinhanca indice->area fica levemente
% menos suave ali. Para ordenar, use sort(caso.catalogo).

caso.catalogo = [ ...
      80.9,   96.6,  108,   121,   130,   153,   156,   178,   182,   199, ...
     254,    325,    332,   414,   419,   453,   523,   533,   705,   733, ...
     809,    862,    989,  1110,  1120,  1230,  1240,  1250,  1270,  1290, ...
    1400,   1440,   1500,  1530,  1550,  1650,  1660,  1780,  1830,  1860, ...
    1910,   2040,   2280,  2230,  2470,  2700,  2760,  3230,  3540,  3600, ...
    4020,   4050,   4280,  4420,  4670,  5360,  5430,  5440,  6380,  7020, ...
    7710,   8040,   8230,  9060,  9380];

assert(numel(caso.catalogo) == 65, ...
    'caso_hadi_20barras:catalogo', ...
    '[HA2003] Sec. 2.2: "In total there are 65 different available sections".');

% -------------------------------------------------------------------------
% SOLUCOES DE REFERENCIA ([HA2003] Tabela 2)
% -------------------------------------------------------------------------
%
% Case 1 = solucao discreta NAO LINEAR   -> 6222 kg
% Case 2 = solucao discreta LINEAR       -> 6253 kg
%
% ================== ATENCAO: DEFEITO DO ARTIGO ==========================
% As areas publicadas na Tabela 2 NAO PERTENCEM ao catalogo de 65 secoes
% publicado na mesma secao. O catalogo termina em 9380 mm^2, mas a tabela
% usa valores maiores:
%
%   Case 1 : A1 = 10400, A3 = 12400, A4 = 10400   (3 dos 10 grupos)
%   Case 2 : A1 = 10400, A3 = 11800, A4 = 10300   (3 dos 10 grupos)
%
% Ou seja, a "solucao discreta" publicada e INATINGIVEL pelo proprio espaco
% de busca que o artigo define. Os outros 7 valores de cada caso estao no
% catalogo normalmente.
%
% CONSEQUENCIA PRATICA: o peso de referencia (6222 / 6253 kg) NAO e um alvo
% valido para este orquestrador, e o desvio contra ele nao mede a qualidade
% da otimizacao. Medido nesta implementacao em 2026-09-02:
%
%   - a solucao publicada do Case 1, avaliada com o solver nao linear deste
%     projeto, da max|sigma| = 160.8 MPa contra o limite de 350 MPa (46% de
%     aproveitamento) e max|u| = 29.1 mm contra 50 mm (58%). Ela deixa folga
%     enorme nas duas restricoes da primeira etapa;
%
%   - um PSO curtissimo (20 particulas, 30 iteracoes) ja encontra projetos
%     VIAVEIS de ~4624 kg usando SOMENTE secoes do catalogo — 26% mais leves
%     que os 6222 kg publicados, com max|sigma| = 265 MPa e max|u| = 48.1 mm,
%     isto e, muito mais proximos do limite.
%
% A leitura mais provavel e que a Tabela 2 reporte a SEGUNDA ETAPA do
% procedimento de [HA2003] (com restricoes de flambagem), e nao a primeira:
% flambagem exige secoes bem maiores nas barras comprimidas, o que explicaria
% ao mesmo tempo as areas grandes e o fato de estourarem o catalogo publicado.
% O artigo nao informa a qual etapa a tabela corresponde.
%
% PORTANTO: compare os resultados deste orquestrador entre si (multi-start,
% sementes, linear vs nao linear), nao contra 6222/6253. As referencias ficam
% registradas para o teste de caracterizacao do peso — que exercita apenas a
% geometria e o agrupamento, e por isso pode usar areas fora do catalogo sem
% problema — e para documentar a discrepancia.
%
% Situacao analoga ja documentada em main_awruch_discreto.m, onde as areas
% de referencia tambem nao pertencem aos catalogos adotados.
% =========================================================================

%                          A1     A2     A3     A4    A5    A6    A7    A8    A9   A10
caso.ref_areas_nao_linear = [10400, 1860, 12400, 10400, 8040, 2700, 1110, 1780, 1290, 8040];
caso.ref_peso_nao_linear  = 6222;    % [kg] Tabela 2, Case 1

caso.ref_areas_linear     = [10400, 2230, 11800, 10300, 8040, 2470,  733, 1650, 6380, 4670];
caso.ref_peso_linear      = 6253;    % [kg] Tabela 2, Case 2

% Aliases genericos (o padrao do orquestrador e a analise nao linear)
caso.ref_areas = caso.ref_areas_nao_linear;
caso.ref_peso  = caso.ref_peso_nao_linear;

% -------------------------------------------------------------------------
% RESTRICAO DE FLAMBAGEM — NAO IMPLEMENTADA, E POR QUE
% -------------------------------------------------------------------------
%
% [HA2003] Sec. 2.2 resolve o problema em duas etapas: a primeira com
% restricoes de tensao e deslocamento, a segunda incluindo flambagem das
% barras. Este orquestrador implementa a PRIMEIRA ETAPA.
%
% A segunda etapa nao e implementavel com os dados disponiveis. As Eqs. (4)
% e (5) de [HA2003] definem a tensao admissivel a compressao em funcao do
% indice de esbeltez S_i = L_i / r_i, onde r_i e o raio de giracao da secao.
% Para os perfis tubulares circulares usados, r depende do diametro E da
% espessura da parede — dois numeros que o artigo nao publica: o catalogo
% da Sec. 2.2 lista somente AREAS. A fonte dessas secoes (Ref. [19] de
% [HA2003]) nao esta no acervo do projeto, em 00_docs/artigos/.
%
% Assumir uma relacao r(A) qualquer produziria numeros com aparencia de
% resultado e sem lastro na fonte. Se a Ref. [19] for obtida, basta
% acrescentar ao catalogo um vetor .raio_giracao e estender avaliar_projeto
% com as Eqs. (3)-(5).

% -------------------------------------------------------------------------
% CONFIGURACAO DAS VARIAVEIS DE PROJETO (para o PSO-RID do Bloco 1)
% -------------------------------------------------------------------------
% 10 variaveis discretas — uma por GRUPO, nao por barra — todas escolhendo
% do mesmo catalogo de 65 secoes.

caso.config_vars = struct('tipo', {}, 'opcoes', {});
for g = 1:numel(caso.grupos)
    caso.config_vars(g).tipo   = 'D';
    caso.config_vars(g).opcoes = caso.catalogo;
end

end


% -------------------------------------------------------------------------
% >>> LOCAL expandir_areas — grupos (10) -> barras (20)
% -------------------------------------------------------------------------

function areas_barras = expandir_areas(areas_grupos, caso)
% EXPANDIR_AREAS  Replica a area de cada grupo nas barras que o compoem.
%
% E a ponte entre o espaco de BUSCA (10 variaveis, o que o PSO enxerga) e o
% espaco de ANALISE (20 barras, o que o FEM exige). Indexacao direta pelo
% mapa inverso montado em caso_hadi_20barras.

areas_barras = areas_grupos(caso.grupo_por_barra);
end


% -------------------------------------------------------------------------
% >>> LOCAL avaliar_projeto — funcao objetivo
% -------------------------------------------------------------------------

function [custo, violacao] = avaliar_projeto(areas_grupos, caso, modo)
% AVALIAR_PROJETO  Peso e violacao de um projeto candidato.
%
% Recebe as 10 areas de GRUPO vindas do PSO, expande para as 20 barras e
% chama o solver do Bloco 2.
%
% A violacao e devolvida SEPARADAMENTE do custo (nao ha penalizacao somada
% ao peso): quem compara e a regra de [DEB2000], que dispensa parametro de
% penalizacao. Ver domina_deb em pso_rid.m.
%
% Medida de violacao: soma linear dos excessos NORMALIZADOS, na forma da
% Eq. (9) de [HA2003]. A justificativa completa esta em avaliar_projeto de
% main_hadi_nao_linear.m.

areas_barras = expandir_areas(areas_grupos, caso);

if strcmp(modo, 'nao_linear')
    [peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas_barras);
else
    [peso, Sigma, u_livre] = fem_linear_solver(caso, areas_barras);
end

violacao = 0;

% [HA2003] Eq. (9) — restricoes normalizadas; ver a justificativa completa
% em avaliar_projeto de main_hadi_nao_linear.m.
g_sigma = abs(Sigma(:))   / caso.sigma_max - 1;
g_desl  = abs(u_livre(:)) / caso.d_max     - 1;

violacao = violacao + sum(max(0, g_sigma));
violacao = violacao + sum(max(0, g_desl));

custo = peso;
end


% -------------------------------------------------------------------------
% >>> LOCAIS de saida
% -------------------------------------------------------------------------

function imprimir_cabecalho(caso, seed, n_runs, pso_params, modo, ref_peso)
n_bits = max(ceil(log2(numel(caso.catalogo))), 1);
fprintf('\n============================================================\n');
fprintf(' %s\n', caso.nome);
fprintf(' Analise: %s\n', upper(strrep(modo, '_', ' ')));
fprintf('============================================================\n');
fprintf(' Barras             : %d, agrupadas em %d variaveis de projeto\n', ...
        size(caso.elements, 2), numel(caso.grupos));
fprintf(' Catalogo           : %d secoes -> %d bits/variavel -> %d dimensoes\n', ...
        numel(caso.catalogo), n_bits, n_bits * numel(caso.grupos));
fprintf(' Espaco de busca    : %d^%d combinacoes\n', numel(caso.catalogo), numel(caso.grupos));
fprintf(' Restricoes         : |sigma| <= %g MPa | |d| <= %g mm\n', ...
        caso.sigma_max, caso.d_max);
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d\n', n_runs);
fprintf(' Particulas         : %d\n', pso_params.n_particulas);
fprintf(' Iteracoes maximas  : %d\n', pso_params.max_iter);
fprintf(' Decodificacao      : %s\n', pso_params.decodificacao_discreta);
fprintf(' Auto-adaptativo    : %d\n', pso_params.auto_adaptativo);
fprintf(' Referencia (artigo): %.2f kg  [NAO e alvo: usa areas fora do catalogo\n', ref_peso);
fprintf('                      e provavelmente inclui flambagem (etapa 2)\n');
fprintf('                      — ver a nota em caso_hadi_20barras]\n');
fprintf('============================================================\n\n');
end


function imprimir_mapa_grupos(areas_grupos, caso)
fprintf(' MAPA DE GRUPOS ([HA2003] Tabela 2)\n');
fprintf(' ------------------------------------------------------------\n');
fprintf(' Grupo |    Area [mm2] | Barras\n');
fprintf(' ------------------------------------------------------------\n');
for g = 1:numel(caso.grupos)
    fprintf('  A%-2d  |  %10.1f | %-18s %s\n', g, areas_grupos(g), ...
            caso.nomes_grupos{g}, mat2str(caso.grupos{g}));
end
fprintf(' ------------------------------------------------------------\n\n');
end


function imprimir_estatisticas_runs(pesos, ref_peso)
validos = pesos(isfinite(pesos));
fprintf(' ESTATISTICAS DAS EXECUCOES\n');
fprintf(' ------------------------------------------------------------\n');
fprintf('  Execucoes         : %d\n', numel(pesos));
fprintf('  Melhor            : %.2f kg\n', min(validos));
fprintf('  Media             : %.2f kg\n', mean(validos));
fprintf('  Desvio padrao     : %.2f kg\n', std(validos));
fprintf('  Pior              : %.2f kg\n', max(validos));
if ~isempty(ref_peso) && ~isnan(ref_peso)
    fprintf('  Gap do melhor     : %+.2f%% em relacao a referencia\n', ...
            100*(min(validos)-ref_peso)/ref_peso);
    fprintf('  (referencia nao e alvo valido — ver nota em caso_hadi_20barras)\n');
end
fprintf(' ------------------------------------------------------------\n\n');
end


function verificar_areas_pertencem_ao_catalogo(areas_grupos, caso)
% Conferencia de coerencia: toda area escolhida deve existir no catalogo.
fora = false;
for g = 1:numel(areas_grupos)
    if ~any(abs(caso.catalogo - areas_grupos(g)) < 1e-9)
        fprintf(2, '  ATENCAO: area do grupo A%d (%.1f mm2) nao pertence ao catalogo.\n', ...
                g, areas_grupos(g));
        fora = true;
    end
end
if ~fora
    fprintf('  OK: todas as %d areas de grupo pertencem ao catalogo de %d secoes.\n\n', ...
            numel(areas_grupos), numel(caso.catalogo));
end
end


function garantir_caminhos()
% Adiciona as pastas do projeto ao path, caso o usuario chame o orquestrador
% diretamente sem ter rodado setup_paths antes.
if exist('pso_rid', 'file') == 2 && exist('fem_nao_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
