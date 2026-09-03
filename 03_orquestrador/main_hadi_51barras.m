function resultado = main_hadi_51barras(seed, n_runs, modo, com_flambagem)
% MAIN_HADI_51BARRAS  Otimizacao da trelica de cobertura de 51 barras, [HA2003] Sec. 6.3.
%
% -------------------------------------------------------------------------
% O QUE ESTE ARQUIVO FAZ
% -------------------------------------------------------------------------
%
% Terceiro exemplo de projeto de [HA2003]. Amarra os mesmos blocos:
%
%   Bloco 1 (otimizador) : pso_rid.m                <- 01_pso_rid/
%   Bloco 2  (solver)    : fem_nao_linear_solver.m  <- 02_fem_nao_linear/
%                          fem_linear_solver.m
%
% O QUE ESTE CASO TRAZ DE NOVO em relacao a main_hadi_20barras:
%
%   1. DUAS HIPOTESES DE CARGA. A estrutura precisa atender tensao e
%      deslocamento sob AMBAS (gravidade e succao de vento). Nenhum caso
%      anterior tinha mais de uma. O solver do Bloco 2 resolve um vetor de
%      cargas por vez, entao avaliar_projeto o chama uma vez por hipotese e
%      soma as violacoes normalizadas — todas as restricoes das duas
%      hipoteses entram na mesma medida, como em [HA2003] Eq. (8)-(9).
%
%   2. OBJETIVO E VOLUME, NAO PESO. A Tabela 3 reporta "Volume (x10^3 mm^3)"
%      e a Sec. 6.3 nao declara densidade nenhuma. Ver a nota sobre .dens
%      em caso_hadi_51barras.
%
%   3. RESTRICAO DE FLAMBAGEM DISPONIVEL. Diferente da Sec. 2.2 (20 barras),
%      aqui [HA2003] fornece o raio de giracao em funcao da area,
%      R = 0.584*A^0.524, o que torna as Eqs. (3)-(5) implementaveis. Por
%      isso este orquestrador cobre os DOIS casos discretos da Tabela 3:
%
%        com_flambagem = false  -> Case 2 (tensao + deslocamento)   PADRAO
%        com_flambagem = true   -> Case 3 (+ flambagem das barras)
%
%   4. Apenas 4 variaveis de projeto para 51 barras (contra 10 para 20).
%
% -------------------------------------------------------------------------
% ESTE ARQUIVO E A FONTE UNICA DO CASO DE 51 BARRAS
% -------------------------------------------------------------------------
%
%   caso = main_hadi_51barras('caso');
%
% Mesmo padrao dos demais orquestradores.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%
%   [HA2003]  Hadi & Alvani (2003), Civil-Comp Press, Paper 37.
%             Sec. 6.3 "Design of 51-bar Roof Truss", Figura 3, Tabela 3.
%   [DF2011]  Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%   [SAKA]    Saka, M.P. — Ref. [20] de [HA2003]. Solucao continua de
%             comparacao e origem da formula R = 0.584*A^0.524.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_hadi_51barras                          % Case 2, nao linear, 5 runs
%   main_hadi_51barras(42, 5, 'linear')         % analise linear
%   main_hadi_51barras(42, 5, 'nao_linear', true)   % Case 3 (com flambagem)
%   caso = main_hadi_51barras('caso');
%
% SAIDA
%   resultado : struct com .melhor_areas (4 grupos) .melhor_areas_barras (51)
%               .melhor_volume .violacao .Sigma .u (por hipotese de carga)
%               .historicos .pesos_por_run .caso .modo .com_flambagem

% -------------------------------------------------------------------------
% ACESSOR DO CASO
% -------------------------------------------------------------------------

if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    resultado = caso_hadi_51barras();
    return;
end

% Acesso as funcoes locais para teste, no mesmo padrao de pso_rid('auxiliares').
% Existe para que 06_testes/ possa verificar a reproducao dos valores
% publicados sem reimplementar a funcao objetivo aqui.
if nargin == 1 && ischar(seed) && strcmp(seed, 'auxiliares')
    resultado = struct( ...
        'expandir_areas',               @expandir_areas,  ...
        'avaliar_projeto',              @avaliar_projeto,  ...
        'resolver_todas_hipoteses',     @resolver_todas_hipoteses,  ...
        'tensao_admissivel_compressao', @tensao_admissivel_compressao);
    return;
end

if nargin < 1 || isempty(seed),          seed          = 42;           end
if nargin < 2 || isempty(n_runs),        n_runs        = 5;            end
if nargin < 3 || isempty(modo),          modo          = 'nao_linear'; end
if nargin < 4 || isempty(com_flambagem), com_flambagem = false;        end

assert(any(strcmp(modo, {'linear', 'nao_linear'})), ...
    'main_hadi_51barras:modoInvalido', ...
    'modo deve ser ''linear'' ou ''nao_linear''. Recebido: ''%s''.', modo);

garantir_caminhos();
rng(seed);

caso = caso_hadi_51barras();

if com_flambagem
    ref_areas  = caso.ref_areas_case3;
    ref_volume = caso.ref_volume_case3;
    rotulo     = 'Case 3 (tensao + deslocamento + flambagem)';
else
    ref_areas  = caso.ref_areas_case2;
    ref_volume = caso.ref_volume_case2;
    rotulo     = 'Case 2 (tensao + deslocamento)';
end

funcao_objetivo = @(a) avaliar_projeto(a, caso, modo, com_flambagem);

% -------------------------------------------------------------------------
% CONFIGURACAO DO PSO-RID
% -------------------------------------------------------------------------
% Espaco de busca pequeno perto dos outros casos: 4 variaveis x 7 bits = 28
% dimensoes, 65^4 = 17.8 milhoes de combinacoes. Enxame e orcamento menores
% bastam, e cada avaliacao custa DOIS solves de FEM (duas hipoteses de carga).

pso_params = struct();
pso_params.n_particulas           = 60;
pso_params.max_iter               = 500;
pso_params.w                      = 1.0;      % [DF2011] Sec. 5
pso_params.c1                     = 1.0;      % [DF2011] Sec. 5
pso_params.c2                     = 2.0;      % [DF2011] Sec. 5
pso_params.pm                     = 0.15;     % [DF2011] Sec. 4.3/5
pso_params.auto_adaptativo        = true;     % [DF2011] Sec. 5
pso_params.decodificacao_discreta = 'proporcional';  % ver [D1] em pso_rid.m
pso_params.tol_estagnacao         = 150;
pso_params.verbose                = true;
pso_params.print_interval         = 50;

imprimir_cabecalho(caso, seed, n_runs, pso_params, modo, rotulo, ref_volume);

% -------------------------------------------------------------------------
% EXECUCAO MULTI-START
% -------------------------------------------------------------------------

historicos      = cell(n_runs, 1);
volumes_por_run = nan(n_runs, 1);
melhor_volume   = inf;
melhor_areas    = [];

for r = 1:n_runs
    fprintf('--- Executando Run %d/%d ---\n', r, n_runs);

    [areas_r, vol_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}      = hist_r;
    volumes_por_run(r) = vol_r;

    if det_r.gbest_viol <= 0 && vol_r < melhor_volume
        melhor_volume = vol_r;
        melhor_areas  = areas_r;
    end

    fprintf('    -> Run %d: volume = %.0f mm3 | violacao = %.3e | iters = %d | avaliacoes = %d\n\n', ...
            r, vol_r, det_r.gbest_viol, det_r.iter_executadas, det_r.n_avaliacoes);
end

if isempty(melhor_areas)
    error('main_hadi_51barras:semSolucaoViavel', ...
        ['Nenhuma das %d execucoes encontrou solucao VIAVEL. ' ...
         'Aumente max_iter/n_particulas ou revise as restricoes.'], n_runs);
end

% -------------------------------------------------------------------------
% AVALIACAO DETALHADA DA MELHOR SOLUCAO
% -------------------------------------------------------------------------

melhor_areas_barras = expandir_areas(melhor_areas, caso);
[vol_final, viol_final] = avaliar_projeto(melhor_areas, caso, modo, com_flambagem);
[Sigma_por_caso, u_por_caso] = resolver_todas_hipoteses(melhor_areas_barras, caso, modo);

imprimir_relatorio(melhor_areas, melhor_areas_barras, vol_final, viol_final, ...
                   ref_areas, ref_volume, Sigma_por_caso, u_por_caso, ...
                   caso, com_flambagem, rotulo);
imprimir_estatisticas_runs(volumes_por_run, ref_volume);
verificar_areas_pertencem_ao_catalogo(melhor_areas, caso);

fig = plot_convergencia(historicos, ref_volume, ...
                        sprintf('PSO-RID — Trelica 51 Barras (%s)', rotulo));
salvar_figura(fig, sprintf('convergencia_hadi_51barras_%s%s', modo, ...
                           repmat('_flambagem', 1, com_flambagem)));

% -------------------------------------------------------------------------
% SAIDA ESTRUTURADA
% -------------------------------------------------------------------------

resultado.caso                = caso;
resultado.modo                = modo;
resultado.com_flambagem       = com_flambagem;
resultado.melhor_areas        = melhor_areas;         % 1 x 4  (grupos)
resultado.melhor_areas_barras = melhor_areas_barras;  % 1 x 51 (barras)
resultado.melhor_volume       = vol_final;
resultado.violacao            = viol_final;
resultado.Sigma               = Sigma_por_caso;       % cell, uma por hipotese
resultado.u                   = u_por_caso;           % cell, uma por hipotese
resultado.historicos          = historicos;
resultado.pesos_por_run       = volumes_por_run;
resultado.seed                = seed;

end


% #########################################################################
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################


% -------------------------------------------------------------------------
% >>> LOCAL caso_hadi_51barras — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_hadi_51barras('caso')
% -------------------------------------------------------------------------

function caso = caso_hadi_51barras()
% CASO_HADI_51BARRAS  Parametros da trelica de cobertura de 51 barras.
%
% -------------------------------------------------------------------------
% REFERENCIA
% -------------------------------------------------------------------------
%
% [HA2003] Sec. 6.3, Figura 3 e Tabela 3 (pagina 13 do PDF em 00_docs/artigos/).
%
% -------------------------------------------------------------------------
% GEOMETRIA (Figura 3)
% -------------------------------------------------------------------------
%
%   Trelica de cobertura de duas aguas, vao 24 m, biapoiada.
%   Banzo inferior HORIZONTAL, nos IMPARES  1,3,5,...,27  (14 nos)
%   Banzo superior INCLINADO,  nos PARES    2,4,6,...,26  (13 nos)
%
%   Malha (Figura 3): 1m + 11 x 2m + 1m = 24 m no banzo inferior; os nos
%   superiores ficam a cada 2 m, deslocados 1 m dos inferiores, o que produz
%   o zigue-zague de 24 diagonais com 1 m de projecao horizontal cada.
%
%                          14 (cume, 2.2 m)
%                        /    \
%             2 __---''''      ''''---__ 26      ^ 1.0 m nas extremidades
%             |                          |       v
%             1---3---5--- ... ---23--25--27
%             ^                          o
%            fixo                      movel
%             |<---------- 24 m ---------->|
%
%   GRUPOS (numeros pequenos na Figura 3):
%     1 = banzo superior (12 barras)      3 = montantes de extremidade (2)
%     2 = banzo inferior (13 barras)      4 = diagonais (24)
%     Total 12 + 13 + 2 + 24 = 51 barras, 4 variaveis de projeto.
%
% -------------------------------------------------------------------------
% ALTURAS: A FIGURA CONTRADIZ O TEXTO, E A FIGURA VENCE
% -------------------------------------------------------------------------
%
% O texto da Sec. 6.3 diz "The slope of the upper chord is taken as 5 [graus]".
% A Figura 3 cota 1.0 m de altura nas extremidades e 1.2 m do cume ate o
% nivel da extremidade — ou seja, 1.2 m de elevacao em 12 m de meio-vao,
% o que da 5.71 graus, nao 5.
%
% As duas hipoteses foram testadas contra os QUATRO volumes publicados na
% Tabela 3 (Volume = soma de A_i * L_i, logo so depende da geometria):
%
%   h_extremidade = 1.0 m, cume = 2.2 m (5.71 graus) -> erro maximo 0.033%
%   h_extremidade = 1.0 m, inclinacao 5.00 graus     -> erro maximo 1.183%
%
% A geometria da figura reproduz os quatro volumes praticamente ao
% arredondamento; a do texto erra 36 vezes mais. Adotada a FIGURA. O "5" do
% texto e provavelmente arredondamento grosseiro de 5.71.
%
% (Um ajuste livre de dois parametros converge para 986 mm e 5.86 graus, com
% erro 0.014% — tao proximo dos valores redondos da figura que nao ha razao
% para preferir numeros quebrados.)
%
% -------------------------------------------------------------------------
% PROPRIEDADES (Sec. 6.3)
% -------------------------------------------------------------------------
%   E = 210x10^3 MPa | deslocamento admissivel 50 mm
%   Tensao de escoamento 350 MPa nos casos DISCRETOS (2 e 3)
%   Catalogo: as mesmas 65 secoes do segundo exemplo (Sec. 2.2)
%
% SAIDA
%   caso : struct com .nodes0 .elements .apoios .F_total .F_casos .E .dens
%          .catalogo .sigma_max .d_max .config_vars .grupos .grupo_por_barra
%          .ref_*

caso.nome = 'Hadi & Alvani (2003) - Trelica de Cobertura 51 Barras (Sec. 6.3)';

% -------------------------------------------------------------------------
% MATERIAL
% -------------------------------------------------------------------------

caso.E = 210e3;          % [MPa]

% DENSIDADE = 1 DE PROPOSITO. O objetivo da Sec. 6.3 e o VOLUME, nao o peso
% (Tabela 3: "Volume (x10^3 mm^3)"), e a secao nao declara densidade alguma.
% Os solvers do Bloco 2 calculam weight = dens * sum(A_i*L_i); com dens = 1
% essa saida passa a ser exatamente o volume em mm^3, sem inventar um valor
% de densidade que o artigo nao fornece.
caso.dens     = 1;
caso.objetivo = 'volume';   % [mm^3]

% -------------------------------------------------------------------------
% RESTRICOES
% -------------------------------------------------------------------------

caso.sigma_max = 350;    % [MPa] escoamento (casos discretos 2 e 3)
caso.d_max     = 50;     % [mm]

% -------------------------------------------------------------------------
% GEOMETRIA
% -------------------------------------------------------------------------

caso.vao        = 24000;   % [mm]
caso.h_apoio    = 1000;    % [mm] altura do banzo superior nas extremidades
caso.h_cume     = 2200;    % [mm] altura do banzo superior no cume

x_inf = [0, 1000, 3000:2000:23000, 24000];          % 14 nos inferiores
x_sup = 0:2000:24000;                               % 13 nos superiores
y_sup = caso.h_apoio + min(x_sup, caso.vao - x_sup) ...
        * (caso.h_cume - caso.h_apoio) / (caso.vao/2);

n_nos       = numel(x_inf) + numel(x_sup);          % 27
caso.nodes0 = zeros(2, n_nos);
caso.nodes0(1, 1:2:end) = x_inf;    % nos IMPARES  = banzo inferior (y = 0)
caso.nodes0(1, 2:2:end) = x_sup;    % nos PARES    = banzo superior
caso.nodes0(2, 2:2:end) = y_sup;

% -------------------------------------------------------------------------
% TOPOLOGIA — 51 barras
% -------------------------------------------------------------------------

banzo_sup = [2:2:24;  4:2:26];    % 12 barras
banzo_inf = [1:2:25;  3:2:27];    % 13 barras
montantes = [1 26;    2 27];      %  2 barras (extremidades)
diagonais = [2:25;    3:26];      % 24 barras (zigue-zague no-a-no)

caso.elements = [banzo_sup, banzo_inf, montantes, diagonais];

assert(size(caso.elements, 2) == 51, ...
    'caso_hadi_51barras:nBarras', '[HA2003] Sec. 6.3: 51 barras.');

% -------------------------------------------------------------------------
% AGRUPAMENTO (numeros pequenos na Figura 3)
% -------------------------------------------------------------------------

n1 = size(banzo_sup, 2); n2 = size(banzo_inf, 2);
n3 = size(montantes, 2); n4 = size(diagonais, 2);

caso.grupos = { ...
    1                     : n1,                ...   % A1 banzo superior
    n1+1                  : n1+n2,             ...   % A2 banzo inferior
    n1+n2+1               : n1+n2+n3,          ...   % A3 montantes extremidade
    n1+n2+n3+1            : n1+n2+n3+n4};            % A4 diagonais

caso.nomes_grupos = {'banzo superior', 'banzo inferior', ...
                     'montantes extr.', 'diagonais'};

caso.grupo_por_barra = zeros(1, size(caso.elements, 2));
for g = 1:numel(caso.grupos)
    caso.grupo_por_barra(caso.grupos{g}) = g;
end
assert(all(caso.grupo_por_barra > 0), ...
    'caso_hadi_51barras:barraSemGrupo', 'Toda barra precisa ter um grupo.');

% -------------------------------------------------------------------------
% CONDICOES DE CONTORNO (Figura 3)
% -------------------------------------------------------------------------
%   no  1 -> apoio FIXO  (GDLs 1 e 2)
%   no 27 -> apoio MOVEL (so vertical: GDL 54; livre em x)
caso.apoios = [1 2 2*27];

% -------------------------------------------------------------------------
% CARREGAMENTO — DUAS HIPOTESES (Figura 3, "Loads in kN")
% -------------------------------------------------------------------------
%
% A Figura 3 desenha o telhado DUAS VEZES, uma para cada hipotese:
%
%   HIPOTESE 1 (gravidade): 10 kN para BAIXO em cada um dos 13 nos do banzo
%   superior, todos iguais, inclusive nas extremidades.
%
%   HIPOTESE 2 (succao/vento): cargas para CIMA, assimetricas —
%     no  2 (extremidade esquerda) ........ 2.5 kN
%     nos 4,6,8,10,12 (agua esquerda) ..... 5 kN cada
%     no 14 (cume) ........................ 5 + 4 = 9 kN
%     nos 16,18,20,22,24 (agua direita) ... 4 kN cada
%     no 26 (extremidade direita) ......... 2 kN
%
% LEITURA DAS SETAS: no cume a figura desenha DUAS setas divergentes (5 e 4)
% porque o no recebe a contribuicao das duas aguas; as demais setas sao
% verticais. Os valores de extremidade (2.5 e 2) sao exatamente METADE de 5 e
% 4, o que confirma area de influencia pela metade nas pontas e fecha a
% leitura: as cargas sao VERTICAIS, nao normais ao plano do telhado (se
% fossem normais, todas as setas apareceriam inclinadas de 5.71 graus, e
% apenas as do cume estao).
%
% As duas hipoteses sao INDEPENDENTES: o projeto precisa atender tensao e
% deslocamento em cada uma delas separadamente (ver avaliar_projeto).

s_dof = 2 * n_nos;
nos_sup = 2:2:26;

F1 = zeros(s_dof, 1);
F1(2*nos_sup) = -10e3;                       % [N] 10 kN para baixo

F2 = zeros(s_dof, 1);
F2(2*2)             =  2.5e3;
F2(2*(4:2:12))      =  5.0e3;
F2(2*14)            =  9.0e3;                % 5 (agua esq.) + 4 (agua dir.)
F2(2*(16:2:24))     =  4.0e3;
F2(2*26)            =  2.0e3;

caso.F_casos     = {F1, F2};
caso.nomes_casos = {'gravidade (10 kN para baixo)', 'succao (5/4 kN para cima)'};
caso.F_total     = F1;   % compatibilidade com quem espera um vetor unico

% -------------------------------------------------------------------------
% CATALOGO — as mesmas 65 secoes da Sec. 2.2
% -------------------------------------------------------------------------
% "Available sections considered for both cases are the same as those used in
%  the second design example above." ([HA2003] Sec. 6.3)
%
% Transcrito da Sec. 2.2, inclusive na ordem publicada (nao monotonica em
% "2280, 2230"). Ver a nota em caso_hadi_20barras.

caso.catalogo = [ ...
      80.9,   96.6,  108,   121,   130,   153,   156,   178,   182,   199, ...
     254,    325,    332,   414,   419,   453,   523,   533,   705,   733, ...
     809,    862,    989,  1110,  1120,  1230,  1240,  1250,  1270,  1290, ...
    1400,   1440,   1500,  1530,  1550,  1650,  1660,  1780,  1830,  1860, ...
    1910,   2040,   2280,  2230,  2470,  2700,  2760,  3230,  3540,  3600, ...
    4020,   4050,   4280,  4420,  4670,  5360,  5430,  5440,  6380,  7020, ...
    7710,   8040,   8230,  9060,  9380];

assert(numel(caso.catalogo) == 65, 'caso_hadi_51barras:catalogo', ...
    'Catalogo de 65 secoes (o mesmo da Sec. 2.2).');

% -------------------------------------------------------------------------
% FLAMBAGEM ([HA2003] Eqs. 3-5) — dados disponiveis NESTE exemplo
% -------------------------------------------------------------------------
%
% "the radius of gyration R is calculated based on the cross-sectional area
%  (A) of the member using R = 0.584*A^0.524, [20]"
%
% Este e o dado que FALTAVA na Sec. 2.2 (20 barras) e que la impediu
% implementar a etapa com flambagem. Aqui ele existe, entao o Case 3 da
% Tabela 3 e reproduzivel. Ver tensao_admissivel_compressao.
%
% RESSALVA: o artigo enuncia a formula ao descrever o Case 1 (secoes de
% cantoneira dupla, continuo). Para os casos discretos as secoes sao os
% perfis tubulares do catalogo da Sec. 2.2, para os quais o artigo nao da
% outra relacao R(A) — e o Case 3 precisa de alguma. Reusa-se a mesma
% formula, que e a unica publicada.

caso.coef_raio_giracao = 0.584;    % R = coef * A^expo
caso.expo_raio_giracao = 0.524;

% -------------------------------------------------------------------------
% SOLUCOES DE REFERENCIA ([HA2003] Tabela 3)
% -------------------------------------------------------------------------
%
%   Saka [20] : solucao CONTINUA de referencia externa
%   Case 1    : continua, flambagem + deslocamento         (nao reproduzivel
%               aqui: variaveis continuas, fora do escopo discreto)
%   Case 2    : DISCRETA, tensao + deslocamento            <- padrao daqui
%   Case 3    : DISCRETA, tensao + deslocamento + flambagem
%
% AO CONTRARIO DA SEC. 2.2, AQUI AS AREAS PUBLICADAS ESTAO NO CATALOGO.
% Case 2 = [1500 1240 523 332] e Case 3 = [1230 1500 1120 533]: todas as oito
% pertencem as 65 secoes. Ou seja, estas referencias SAO alvos legitimos —
% diferente do caso de 20 barras, onde a Tabela 2 usava areas inexistentes.

caso.ref_areas_saka   = [1673, 1149,  378, 721];   caso.ref_volume_saka   = 101537e3;
caso.ref_areas_case1  = [1473, 1095, 1550, 665];   caso.ref_volume_case1  =  95172e3;
caso.ref_areas_case2  = [1500, 1240,  523, 332];   caso.ref_volume_case2  =  82098e3;
caso.ref_areas_case3  = [1230, 1500, 1120, 533];   caso.ref_volume_case3  =  92169e3;

% Aliases genericos (o padrao do orquestrador e o Case 2)
caso.ref_areas = caso.ref_areas_case2;
caso.ref_peso  = caso.ref_volume_case2;

% -------------------------------------------------------------------------
% DISCREPANCIA DOCUMENTADA: OS OTIMOS PUBLICADOS FICAM ~7% ACIMA DO LIMITE
% DE DESLOCAMENTO NESTE MODELO   (investigacao de 2026-09-02)
% -------------------------------------------------------------------------
%
% Avaliando as areas publicadas com o FEM deste projeto, sob a hipotese de
% carga gravitacional:
%
%   Case 2 [1500 1240 523 332] : max|sigma| = 213.0 MPa (61% de 350)  OK
%                                max|u|     =  53.53 mm (107% de 50)  ESTOURA
%   Case 3 [1230 1500 1120 533]: max|u|     =  50.97 mm (102% de 50)  ESTOURA
%
% A TENSAO SOBRA nos dois; o que estoura e SO o deslocamento, e por pouco.
% Note que a restricao esta ATIVA nos dois projetos (53.5 e 51.0 contra 50),
% que e exatamente o que se espera de um otimo — sinal de que o modelo esta
% essencialmente certo e a diferenca e de poucos por cento de rigidez, nao
% erro de topologia.
%
% O QUE FOI DESCARTADO como causa:
%
%   - CARGA NAS EXTREMIDADES. Reduzir para 5 kN (metade da area de
%     influencia) nos nos 2 e 26 muda max|u| de 53.53 para 53.44 mm: os nos
%     de extremidade ficam sobre os apoios e quase nao fletem. Nao explica.
%
%   - APOIOS. Com dois apoios FIXOS o deslocamento cai para 31.85 mm, mas ai
%     os otimos publicados ficariam a 64% do limite — longe da fronteira, o
%     que nao combina com solucao otima. Alem disso a Figura 3 desenha
%     claramente rolete no no 27. Descartado.
%
%   - GEOMETRIA. Varredura de h_apoio e h_cume: so ha viabilidade a partir de
%     h_cume >= 2300 mm, e isso degrada o erro nos quatro volumes da Tabela 3
%     de 0.033% para 0.77%..1.7% — 20 a 50 vezes pior. O compromisso e
%     monotono: toda geometria que viabiliza os otimos piora o volume. Como o
%     volume e a UNICA grandeza da Tabela 3 que da para conferir de forma
%     independente, e ele confirma 1000/2200 mm com folga, a geometria fica.
%
% CONCLUSAO: mantem-se a leitura fiel da figura, e a diferenca de ~7% na
% rigidez fica registrada como divergencia em aberto contra [HA2003]. Ela nao
% invalida o caso — o otimizador daqui encontra solucoes VIAVEIS com volume
% 6% a 8% acima do publicado, coerente com um modelo um pouco mais flexivel
% que o do artigo. Compare execucoes deste orquestrador entre si; use a
% Tabela 3 como ordem de grandeza, nao como alvo exato.
%
% (O que NAO foi possivel testar: o artigo nao publica deslocamentos nem
% tensoes das solucoes da Tabela 3, so areas e volume. Sem isso nao da para
% isolar se a diferenca esta na analise, na medida do deslocamento ou no
% proprio resultado publicado.)

% -------------------------------------------------------------------------
% VARIAVEIS DE PROJETO — 4 discretas, uma por grupo
% -------------------------------------------------------------------------

caso.config_vars = struct('tipo', {}, 'opcoes', {});
for g = 1:numel(caso.grupos)
    caso.config_vars(g).tipo   = 'D';
    caso.config_vars(g).opcoes = caso.catalogo;
end

end


% -------------------------------------------------------------------------
% >>> LOCAL expandir_areas — grupos (4) -> barras (51)
% -------------------------------------------------------------------------

function areas_barras = expandir_areas(areas_grupos, caso)
% Ponte entre o espaco de BUSCA (4 variaveis) e o de ANALISE (51 barras).
areas_barras = areas_grupos(caso.grupo_por_barra);
end


% -------------------------------------------------------------------------
% >>> LOCAL resolver_todas_hipoteses
% -------------------------------------------------------------------------

function [Sigma_por_caso, u_por_caso, L0] = resolver_todas_hipoteses(areas_barras, caso, modo)
% Roda o solver do Bloco 2 uma vez por hipotese de carga.
%
% O solver le o carregamento de problema.F_total, entao trocamos esse campo a
% cada hipotese. Tudo o mais (geometria, apoios, material) e identico.

n_casos        = numel(caso.F_casos);
Sigma_por_caso = cell(n_casos, 1);
u_por_caso     = cell(n_casos, 1);
prob           = caso;

for k = 1:n_casos
    prob.F_total = caso.F_casos{k};
    if strcmp(modo, 'nao_linear')
        [~, S, u, d] = fem_nao_linear_solver(prob, areas_barras);
    else
        [~, S, u, d] = fem_linear_solver(prob, areas_barras);
    end
    Sigma_por_caso{k} = S;
    u_por_caso{k}     = u;
    L0                = d.L0;
end
end


% -------------------------------------------------------------------------
% >>> LOCAL tensao_admissivel_compressao — [HA2003] Eqs. (4)-(5)
% -------------------------------------------------------------------------

function sigma_a = tensao_admissivel_compressao(areas_barras, L0, caso)
% Tensao admissivel a COMPRESSAO por barra, pela especificacao AISC citada em
% [HA2003] Eqs. (4) e (5).
%
%   R = 0.584 * A^0.524                 raio de giracao ([HA2003] Sec. 6.3)
%   S = L / R                           indice de esbeltez
%   C = sqrt(2*pi^2*E / Fy)
%
%   S <  C :  sigma_a = [1 - S^2/(2C^2)] * Fy / [5/3 + 3S/(8C) - S^3/(8C^3)]
%   S >= C :  sigma_a = 12*pi^2*E / (23*S^2)

Fy = caso.sigma_max;
E  = caso.E;

R = caso.coef_raio_giracao * areas_barras(:).^caso.expo_raio_giracao;
S = L0(:) ./ R;
C = sqrt(2 * pi^2 * E / Fy);

sigma_a = zeros(size(S));

curto = S < C;
sigma_a(curto) = (1 - S(curto).^2 ./ (2*C^2)) * Fy ./ ...
                 (5/3 + 3*S(curto)./(8*C) - S(curto).^3 ./ (8*C^3));

longo = ~curto;
sigma_a(longo) = 12 * pi^2 * E ./ (23 * S(longo).^2);
end


% -------------------------------------------------------------------------
% >>> LOCAL avaliar_projeto — funcao objetivo
% -------------------------------------------------------------------------

function [custo, violacao] = avaliar_projeto(areas_grupos, caso, modo, com_flambagem)
% AVALIAR_PROJETO  Volume e violacao de um projeto candidato.
%
% CUSTO = VOLUME em mm^3 (objetivo da Sec. 6.3). Como caso.dens = 1, a saida
% "weight" dos solvers do Bloco 2 ja e o volume — ver a nota em
% caso_hadi_51barras.
%
% VIOLACAO: soma dos excessos NORMALIZADOS ([HA2003] Eq. 9, ver [D8] no
% relatorio de decisoes), acumulada sobre AS DUAS HIPOTESES DE CARGA. Cada
% hipotese contribui com suas n tensoes e m deslocamentos; o projeto so e
% viavel se atender as duas. Devolvida SEPARADA do custo — [DEB2000].
%
% Com com_flambagem = true, a tensao admissivel deixa de ser Fy nos dois
% sentidos e passa a ser ([HA2003] Eq. 3):
%     tracao     : 0.6 * Fy                       (especificacao AISC)
%     compressao : Eqs. (4)-(5), funcao da esbeltez de cada barra

areas_barras = expandir_areas(areas_grupos, caso);

prob     = caso;
volume   = NaN;
violacao = 0;

for k = 1:numel(caso.F_casos)
    prob.F_total = caso.F_casos{k};

    if strcmp(modo, 'nao_linear')
        [volume, Sigma, u_livre, diag] = fem_nao_linear_solver(prob, areas_barras);
    else
        [volume, Sigma, u_livre, diag] = fem_linear_solver(prob, areas_barras);
    end

    % --- restricao de tensao ---
    if com_flambagem
        % [HA2003] Eq. (3): admissiveis distintos em tracao e compressao
        sig_comp = tensao_admissivel_compressao(areas_barras, diag.L0, caso);
        sig_trac = 0.6 * caso.sigma_max;

        limite = sig_trac * ones(size(Sigma(:)));
        limite(Sigma(:) < 0) = sig_comp(Sigma(:) < 0);
    else
        % Case 2: escoamento nos dois sentidos ([HA2003] Eq. 2)
        limite = caso.sigma_max * ones(size(Sigma(:)));
    end

    g_sigma = abs(Sigma(:)) ./ limite - 1;
    g_desl  = abs(u_livre(:)) / caso.d_max - 1;

    violacao = violacao + sum(max(0, g_sigma));
    violacao = violacao + sum(max(0, g_desl));
end

custo = volume;   % identico nas duas hipoteses (so depende de A e L0)
end


% -------------------------------------------------------------------------
% >>> LOCAIS de saida
% -------------------------------------------------------------------------

function imprimir_cabecalho(caso, seed, n_runs, pso_params, modo, rotulo, ref_volume)
n_bits = max(ceil(log2(numel(caso.catalogo))), 1);
fprintf('\n============================================================\n');
fprintf(' %s\n', caso.nome);
fprintf(' %s | Analise: %s\n', rotulo, upper(strrep(modo, '_', ' ')));
fprintf('============================================================\n');
fprintf(' Barras             : %d, agrupadas em %d variaveis de projeto\n', ...
        size(caso.elements, 2), numel(caso.grupos));
fprintf(' Hipoteses de carga : %d (avaliadas a cada chamada da f. objetivo)\n', ...
        numel(caso.F_casos));
for k = 1:numel(caso.F_casos)
    fprintf('    %d) %s\n', k, caso.nomes_casos{k});
end
fprintf(' Catalogo           : %d secoes -> %d bits/variavel -> %d dimensoes\n', ...
        numel(caso.catalogo), n_bits, n_bits * numel(caso.grupos));
fprintf(' Objetivo           : VOLUME [mm3]  (Tabela 3 de [HA2003])\n');
fprintf(' Restricoes         : |d| <= %g mm | tensao: %s\n', caso.d_max, ...
        ternario(strcmp(rotulo(1:6),'Case 3'), ...
                 'AISC Eqs. (3)-(5), com flambagem', ...
                 sprintf('|sigma| <= %g MPa', caso.sigma_max)));
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d\n', n_runs);
fprintf(' Particulas         : %d | Iteracoes maximas: %d\n', ...
        pso_params.n_particulas, pso_params.max_iter);
fprintf(' Referencia (artigo): %.0f mm3 (%.0f x10^3)\n', ref_volume, ref_volume/1e3);
fprintf('============================================================\n\n');
end


function s = ternario(cond, a, b)
if cond, s = a; else, s = b; end
end


function imprimir_relatorio(areas_g, areas_b, volume, viol, ref_areas, ref_volume, ...
                            Sigma_por_caso, u_por_caso, caso, com_flambagem, rotulo)
fprintf('\n============================================================\n');
fprintf(' RELATORIO: TRELICA 51 BARRAS — %s\n', upper(rotulo));
fprintf('============================================================\n');
fprintf(' AREAS POR GRUPO [mm2]\n');
fprintf(' ------------------------------------------------------------\n');
fprintf(' Grupo | %-17s |   PSO-RID | Referencia | Barras\n', 'Descricao');
fprintf(' ------------------------------------------------------------\n');
for g = 1:numel(areas_g)
    fprintf('  A%-2d  | %-17s | %9.1f | %10.1f | %d\n', g, caso.nomes_grupos{g}, ...
            areas_g(g), ref_areas(g), numel(caso.grupos{g}));
end
fprintf(' ------------------------------------------------------------\n');
fprintf(' Volume PSO-RID     : %.0f mm3  (%.0f x10^3)\n', volume, volume/1e3);
fprintf(' Volume referencia  : %.0f mm3  (%.0f x10^3)\n', ref_volume, ref_volume/1e3);
fprintf(' Diferenca          : %+.2f%%\n', 100*(volume-ref_volume)/ref_volume);
fprintf(' Violacao total     : %.6e\n', viol);
fprintf(' ------------------------------------------------------------\n');

if com_flambagem
    sig_comp = tensao_admissivel_compressao(areas_b, ...
                   comprimentos_barras(caso), caso);
end

for k = 1:numel(Sigma_por_caso)
    S = Sigma_por_caso{k};  u = u_por_caso{k};
    fprintf(' HIPOTESE %d — %s\n', k, caso.nomes_casos{k});
    if com_flambagem
        lim = 0.6*caso.sigma_max * ones(size(S(:)));
        lim(S(:) < 0) = sig_comp(S(:) < 0);
        pior = max(abs(S(:)) ./ lim);
        fprintf('   tensao   : pior aproveitamento %.1f%% do admissivel (variavel por barra)\n', 100*pior);
    else
        fprintf('   tensao   : max |sigma| = %8.2f MPa  (limite %g)  -> %.1f%%\n', ...
                max(abs(S)), caso.sigma_max, 100*max(abs(S))/caso.sigma_max);
    end
    fprintf('   desloc.  : max |d|     = %8.2f mm   (limite %g)  -> %.1f%%\n', ...
            max(abs(u)), caso.d_max, 100*max(abs(u))/caso.d_max);
end
fprintf(' ------------------------------------------------------------\n\n');
end


function L = comprimentos_barras(caso)
L = zeros(size(caso.elements, 2), 1);
for e = 1:numel(L)
    L(e) = norm(caso.nodes0(:, caso.elements(2,e)) - caso.nodes0(:, caso.elements(1,e)));
end
end


function imprimir_estatisticas_runs(volumes, ref_volume)
validos = volumes(isfinite(volumes));
fprintf(' ESTATISTICAS DAS EXECUCOES  (volume em mm3)\n');
fprintf(' ------------------------------------------------------------\n');
fprintf('  Execucoes         : %d\n', numel(volumes));
fprintf('  Melhor            : %.0f\n', min(validos));
fprintf('  Media             : %.0f\n', mean(validos));
fprintf('  Desvio padrao     : %.0f\n', std(validos));
fprintf('  Pior              : %.0f\n', max(validos));
if ~isempty(ref_volume) && ~isnan(ref_volume)
    fprintf('  Gap do melhor     : %+.2f%% em relacao a referencia\n', ...
            100*(min(validos)-ref_volume)/ref_volume);
end
fprintf(' ------------------------------------------------------------\n\n');
end


function verificar_areas_pertencem_ao_catalogo(areas_grupos, caso)
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
if exist('pso_rid', 'file') == 2 && exist('fem_nao_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
