function resultado = main_hadi_nao_linear(seed, n_runs, n_workers)
arguments
    % seed aceita tambem char/string por causa da chamada-sentinela
    % main_hadi_nao_linear('caso') — ver mais abaixo.
    seed {mustBeA(seed, ["double","char","string"])} = 42
    n_runs    (1,1) double {mustBePositive, mustBeInteger} = 5
    n_workers (1,1) double {mustBeNonnegative, mustBeInteger} = 0
end
% MAIN_HADI_NAO_LINEAR  Otimizacao da trelica de 10 barras com analise NAO LINEAR.
%
% -------------------------------------------------------------------------
% O QUE ESTE ARQUIVO FAZ
% -------------------------------------------------------------------------
%
% Este e o unico arquivo que precisa ser lido para entender o experimento
% completo. Ele amarra os blocos:
%
%   Bloco 1 (otimizador) : pso_rid.m                <- 01_pso_rid/
%   Bloco 2  (solver)    : fem_nao_linear_solver.m  <- 02_fem_nao_linear/
%
% Toda a logica especifica do experimento (DEFINICAO DO PROBLEMA, funcao
% objetivo, laco multi-start, selecao da melhor solucao) esta embutida aqui
% como funcoes locais no fim do arquivo. Ficam fora apenas os blocos longos
% de impressao e grafico, que poluiriam a leitura do fluxo principal:
%
%   03_orquestrador/auxiliares/relatorio_comparativo.m
%   03_orquestrador/auxiliares/plot_convergencia.m
%   03_orquestrador/auxiliares/salvar_figura.m
%
% FONTE UNICA DO CASO HADI: geometria, material, cargas, catalogo e solucao
% de referencia vivem na funcao local caso_hadi_10barras, exposta por
% main_hadi_nao_linear('caso'). Usam esse acessor main_hadi_linear.m,
% main_estudo_estatistico.m, main_awruch_discreto.m (herda a geometria) e os
% testes do Bloco 2.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%
%   [HA2003]  Hadi & Alvani (2003), Civil-Comp Press, Paper 37.
%             Benchmark de 10 barras, Sec. 6.1 e Tabela 1 (Case 2, discreto).
%   [DF2011]  Datta & Figueira, Applied Soft Computing 11 (2011) 3625-3633.
%             Formulacao do PSO real-inteiro-discreto.
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%             Tratamento de restricoes sem parametro de penalizacao.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_hadi_nao_linear              % semente 42, 5 execucoes
%   main_hadi_nao_linear(123)         % outra semente
%   main_hadi_nao_linear(42, 10)      % 10 execucoes
%   main_hadi_nao_linear(42, 5, 3)    % 10 execucoes, ate 3 workers em paralelo
%   r = main_hadi_nao_linear;         % devolve struct com os resultados
%
% n_workers (padrao 0 = serial) paraleliza o laco multi-start via parfor.
% Efeito colateral aceito: cada run passa a ter semente PROPRIA
% (seed + run - 1) em vez de todas consumirem o mesmo stream em sequencia —
% necessario para o parfor ser seguro, e documentado na Secao 4 abaixo.
% preparar_pool.m aplica um teto por MEMORIA disponivel, nao por nucleos.
%
% SAIDA
%   resultado : struct com .melhor_areas .melhor_peso .violacao .Sigma .u
%               .historicos .pesos_por_run .caso
%
% See also pso_rid, fem_nao_linear_solver, main_hadi_linear

% -------------------------------------------------------------------------
% ACESSOR DO CASO: main_hadi_nao_linear('caso') devolve so o struct do
% problema, sem rodar nenhuma otimizacao (ver cabecalho).
% -------------------------------------------------------------------------

if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    resultado = caso_hadi_10barras();
    return;
end

garantir_caminhos({'pso_rid', 'fem_nao_linear_solver'});
n_workers = preparar_pool(n_workers, n_runs);

% -------------------------------------------------------------------------
% 1. PROBLEMA — apenas parametros, nenhuma logica (funcao local no fim)
% -------------------------------------------------------------------------

caso = caso_hadi_10barras();

% -------------------------------------------------------------------------
% 2. FUNCAO OBJETIVO
% -------------------------------------------------------------------------
% Acopla o solver FEM (Bloco 2a) as restricoes do problema. Definida como
% funcao local no fim deste arquivo (avaliar_projeto)
% -------------------------------------------------------------------------

funcao_objetivo = @(areas) avaliar_projeto(areas, caso, 'nao_linear');

% -------------------------------------------------------------------------
% 3. CONFIGURACAO DO PSO-RID (Bloco 1)
% -------------------------------------------------------------------------
% Valores iniciais conforme [DF2011] Sec. 5, pag. 3628. Com
% auto_adaptativo = true, cada parametro e sorteado a cada iteracao no
% intervalo [0, valor_inicial].
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
pso_params.tol_estagnacao         = 200;
pso_params.verbose                = true;
pso_params.print_interval         = 100;

% -------------------------------------------------------------------------
% 4. EXECUCAO MULTI-START
% -------------------------------------------------------------------------

imprimir_cabecalho(caso, seed, n_runs, pso_params, 'NAO LINEAR');

historicos      = cell(n_runs, 1);
pesos_por_run   = nan(n_runs, 1);
areas_por_run   = cell(n_runs, 1);
viaveis_por_run = false(n_runs, 1);
linhas_log      = cell(n_runs, 1);

% Cada run e independente. Rodar em parfor (n_workers=0 -> serial, o padrao)
% exige que cada run tenha sua PROPRIA semente explicita, em vez de todas
% consumirem em sequencia o mesmo stream global — o mesmo motivo, e a mesma
% solucao, ja documentados em pso_rid.m (Secao C) e no cabecalho deste
% arquivo. Consequencia: o resultado da run r passa a depender so de "seed"
% e "r", nao de quantas runs vieram antes.
%
% [NOTA] fprintf de dentro de pso_rid (pso_params.verbose) ainda pode
% intercalar entre workers quando n_workers>0 — e cosmetico, nao afeta o
% resultado. As linhas de resumo abaixo (linhas_log) sao guardadas e
% impressas em ordem DEPOIS do laco, como em main_estudo_estatistico.m.

parfor (r = 1:n_runs, n_workers)
    rng(seed + r - 1, 'twister');   % semente propria por run — parfor-safe

    [areas_r, peso_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}      = hist_r;
    pesos_por_run(r)   = peso_r;
    areas_por_run{r}   = areas_r;
    viaveis_por_run(r) = (det_r.gbest_viol <= 0);

    linhas_log{r} = sprintf(['--- Run %d/%d ---\n    -> peso = %.2f kg | ' ...
        'violacao = %.3e | iters = %d | avaliacoes FEM = %d\n\n'], ...
        r, n_runs, peso_r, det_r.gbest_viol, det_r.iter_executadas, det_r.n_avaliacoes);
end

fprintf('%s', linhas_log{:});

% Seleciona o melhor apenas entre solucoes VIAVEIS ([DEB2000] criterio 1).
% Fica em laco serial separado — nao pode entrar no parfor porque tem
% dependencia de ORDEM entre iteracoes (reducao sequencial).
melhor_peso  = inf;
melhor_areas = [];
for r = 1:n_runs
    if viaveis_por_run(r) && pesos_por_run(r) < melhor_peso
        melhor_peso  = pesos_por_run(r);
        melhor_areas = areas_por_run{r};
    end
end

if isempty(melhor_areas)
    error('main_hadi_nao_linear:semSolucaoViavel', ...
        ['Nenhuma das %d execucoes encontrou solucao VIAVEL. ' ...
         'Aumente max_iter/n_particulas ou revise as restricoes.'], n_runs);
end

% -------------------------------------------------------------------------
% 5. AVALIACAO DETALHADA DA MELHOR SOLUCAO
% -------------------------------------------------------------------------

[peso_final, Sigma_final, u_final] = fem_nao_linear_solver(caso, melhor_areas);
[~, viol_final] = avaliar_projeto(melhor_areas, caso, 'nao_linear');

% -------------------------------------------------------------------------
% 6. RELATORIO E GRAFICOS (auxiliares — blocos longos de saida)
% -------------------------------------------------------------------------

relatorio_comparativo('Trelica 10 Barras Hadi (Nao Linear)', ...
                      melhor_areas, peso_final, ...
                      caso.ref_areas, caso.ref_peso, ...
                      Sigma_final, caso.sigma_max, ...
                      u_final, caso.d_max, viol_final);

imprimir_estatisticas_runs(pesos_por_run, caso.ref_peso);

fig = plot_convergencia(historicos, caso.ref_peso, ...
                        'PSO-RID — Trelica 10 Barras (Analise Nao Linear)');
salvar_figura(fig, 'convergencia_hadi_nao_linear');

% -------------------------------------------------------------------------
% 7. SAIDA ESTRUTURADA
% -------------------------------------------------------------------------

resultado.caso          = caso;
resultado.melhor_areas  = melhor_areas;
resultado.melhor_peso   = peso_final;
resultado.violacao      = viol_final;
resultado.Sigma         = Sigma_final;
resultado.u             = u_final;
resultado.historicos    = historicos;
resultado.pesos_por_run = pesos_por_run;
resultado.seed          = seed;

end


% #########################################################################
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################


% -------------------------------------------------------------------------
% >>> LOCAL caso_hadi_10barras — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_hadi_nao_linear('caso')
% -------------------------------------------------------------------------

function caso = caso_hadi_10barras()
% CASO_HADI_10BARRAS  Parametros do benchmark de trelica plana de 10 barras.
%
% Sem logica de solver; ver a nota no cabecalho deste arquivo.
%
% -------------------------------------------------------------------------
% REFERENCIA
% -------------------------------------------------------------------------
%
% [HA2003] Hadi, M.N.S.; Alvani, K.S. "Discrete Optimum Design of
%          Geometrically Non-Linear Trusses using Genetic Algorithms".
%          Civil-Comp Press, 2003, Paper 37.  ->  Secao 6.1, Figura 1, Tabela 1
%          -> 00_docs/artigos/Discrete Optimum Design of Geometrically
%             Non-Linear Trusses using Genetic Algorithms.PDF
%
% -------------------------------------------------------------------------
% GEOMETRIA (Figura 1 de [HA2003])
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% PROPRIEDADES (Secao 6.1 de [HA2003])
% -------------------------------------------------------------------------
%   "the modulus of elasticity was 6.89 x 10^4 MPa, and the material density
%    2770 kg/m^3 the allowable displacement was limited to 50.8 mm and the
%    allowable stress to 172.25 MPa in both tension and compression for all
%    members. The problem has 10 design variables."
%
% UNIDADES: mm / N / MPa / kg  (densidade convertida para kg/mm^3)
%
% SAIDA
%   caso : struct com os campos exigidos pelos solvers do Bloco 2
%          (.nodes0 .elements .apoios .F_total .E .dens) mais os dados de
%          otimizacao (.catalogo .sigma_max .d_max .config_vars .ref_*)

caso.nome = 'Hadi & Alvani (2003) - Trelica 10 Barras (Case 2, discreto)';

% -------------------------------------------------------------------------
% O PROBLEMA E ORIGINALMENTE IMPERIAL  (investigacao de 2026-08-31)
% -------------------------------------------------------------------------
%
% [HA2003] enuncia tudo em SI, mas todos os valores sao conversoes de
% numeros redondos em unidades imperiais:
%
%   9144 mm      = 360 in          (modulo da malha)
%   50.80 mm     = 2 in            (deslocamento admissivel)
%   172.25 MPa   = 25 ksi          (tensao admissivel)
%   6.89e4 MPa   = 10^4 ksi        (arredondado de 68947.573 MPa)
%   2770 kg/m^3  = 0.1 lb/in^3     (arredondado de 2767.990 kg/m^3)
%   445.4 kN     = 100 kip         (arredondado de 444.822 kN)
%   catalogo     = 0.1, 1, 2, 5, 8, 12, 15, 18, 20, 25, 30, 35, 45 in^2
%
% ISSO EXPLICA A DISCREPANCIA DE PESO que ficou em aberto desde 2026-08-28:
% o peso publicado (2325.2 kg) foi calculado com rho = 0.1 lb/in^3 EXATO,
% nao com o 2770 kg/m^3 que o proprio texto do artigo declara. Usando a
% densidade exata, esta implementacao reproduz a Tabela 1 assim:
%
%   dens = 2770    (texto do artigo) -> 2326.942 kg   (+0.075% vs publicado)
%   dens = 2767.99 (0.1 lb/in^3)     -> 2325.254 kg   (+0.002% vs publicado)
%
% Verificado contra as SEIS colunas da Tabela 1 com areas listadas: cinco
% reproduzem o peso publicado com erro <= 0.002%. Ver o relatorio em
% 00_docs/notas_e_relatorios/ e o teste de caracterizacao em
% 06_testes/TestFemNaoLinear.m.
%
% ESCOLHA DE PADRAO: mantem-se os valores DECLARADOS NO TEXTO do artigo
% (2770 e 6.89e4), para que os resultados ja registrados em 04_resultados/
% continuem reproduziveis. Os valores exatos ficam disponiveis nos campos
% .dens_exata e .E_exato para quem quiser reproduzir a Tabela 1 ao milesimo.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% MATERIAL
% -------------------------------------------------------------------------

caso.E    = 6.89e4;      % [MPa]      modulo de elasticidade  [HA2003] Sec. 6.1
caso.dens = 2.770e-6;    % [kg/mm^3]  2770 kg/m^3 = 2.770e-6 kg/mm^3

% Equivalentes imperiais exatos (ver nota acima). Fatores de conversao
% exatos por definicao: 1 lb = 0.45359237 kg, 1 in = 25.4 mm.
caso.E_exato     = 1e4 * 6.894757293168361;              % [MPa]     10^4 ksi
caso.dens_exata  = 0.1 * 0.45359237 / 16.387064e-6 * 1e-9;  % [kg/mm^3] 0.1 lb/in^3

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

%   barra        1  2  3  4  5  6  7  8  9 10
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
%
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


function [custo, violacao] = avaliar_projeto(areas, caso, tipo_analise)
% AVALIAR_PROJETO  Funcao objetivo: peso sujeito a tensao e deslocamento.
%
% Formulacao do problema — [HA2003] Eq. (1) e (2):
%
%   Minimizar   W(A) = sum_i  A_i * L_i * rho                        Eq. (1)
%   sujeito a   |sigma_i| <= sigma_admissivel      i = 1..n          Eq. (2)
%               |d_j|     <= d_admissivel          j = 1..m
%
% A violacao e devolvida SEPARADAMENTE do custo (nao ha penalizacao somada
% ao peso), porque o PSO-RID trata restricoes pela regra de [DEB2000], que
% dispensa parametro de penalizacao. Ver a funcao local domina_deb,
% no fim de 01_pso_rid/pso_rid.m.
%
% Medida de violacao: soma linear dos excessos NORMALIZADOS, na forma da
% Eq. (9) de [HA2003]. A justificativa esta no corpo da funcao.

if strcmp(tipo_analise, 'linear')
    [peso, Sigma, u_livre] = fem_linear_solver(caso, areas);
else
    [peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);
end

violacao = 0;

% -------------------------------------------------------------------------
% [HA2003] Eq. (9): restricoes na forma NORMALIZADA e adimensional
%
%   g_sigma_i = |sigma_i| / sigma_adm - 1 <= 0        i = 1..n_barras
%   g_d_j     = |d_j|     / d_adm     - 1 <= 0        j = 1..n_gdl_livres
%
% Violacao total = soma dos excessos, <g> = max(0, g) ([DEB2000] pag. 316).
%
% POR QUE NORMALIZAR: sem dividir pelo admissivel, os dois excessos entram
% na soma em unidades diferentes — MPa para tensao e mm para deslocamento —
% e a razao entre eles passa a depender da escolha de unidades, nao da
% fisica. Com sigma_adm = 172.25 MPa e d_adm = 50.80 mm, 1 mm de excesso de
% deslocamento pesava o mesmo que 1 MPa de excesso de tensao, embora o
% primeiro consuma 2.0% da folga e o segundo 0.6%. Normalizado, "1.0" quer
% dizer a mesma coisa nas duas familias: 100% de excesso sobre o limite.
% Isso importa porque o criterio 3 de [DEB2000] ordena inviaveis SOMENTE
% pela violacao — uma medida incomensuravel ordena pela unidade escolhida.
%
% POR QUE LINEAR E NAO QUADRATICO: a versao anterior somava os quadrados
% dos excessos brutos. O quadrado amplifica a incomensurabilidade em vez de
% corrigi-la (MPa^2 contra mm^2) e nao consta de [HA2003] nem de [DEB2000].
% -------------------------------------------------------------------------

g_sigma = abs(Sigma(:))   / caso.sigma_max - 1;
g_desl  = abs(u_livre(:)) / caso.d_max     - 1;

violacao = violacao + sum(max(0, g_sigma));
violacao = violacao + sum(max(0, g_desl));

custo = peso;
end


function imprimir_cabecalho(caso, seed, n_runs, pso_params, tipo)
fprintf('\n============================================================\n');
fprintf(' OTIMIZACAO PSO-RID — %s\n', caso.nome);
fprintf(' Analise estrutural: %s\n', tipo);
fprintf('------------------------------------------------------------\n');
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d\n', n_runs);
fprintf(' Particulas         : %d\n', pso_params.n_particulas);
fprintf(' Iteracoes maximas  : %d\n', pso_params.max_iter);
fprintf(' Decodificacao      : %s\n', pso_params.decodificacao_discreta);
fprintf(' Auto-adaptativo    : %d\n', pso_params.auto_adaptativo);
fprintf(' Referencia (artigo): %.2f kg\n', caso.ref_peso);
fprintf('============================================================\n\n');
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
end
fprintf(' ------------------------------------------------------------\n\n');
end


% garantir_caminhos e preparar_pool sao helpers compartilhados em
% 03_orquestrador/auxiliares/ (ver garantir_caminhos.m e preparar_pool.m) —
% nao ha mais funcao local aqui.
