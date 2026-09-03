function resultado = main_datta_mola(seed, n_runs)
% MAIN_DATTA_MOLA  Projeto da mola helicoidal de compressao, [DF2011] Sec. 5.2.
%
% -------------------------------------------------------------------------
% O QUE ESTE ARQUIVO FAZ
% -------------------------------------------------------------------------
%
% Segundo problema analitico do projeto (nao usa o Bloco 2 / FEM). Tudo mora
% neste arquivo: definicao do problema, funcao objetivo, restricoes e laco de
% execucao, como funcoes locais no fim.
%
%   Bloco 1 (otimizador) : pso_rid.m   <- 01_pso_rid/
%   Resto                : local, aqui
%
% -------------------------------------------------------------------------
% POR QUE ESTE E O CASO MAIS IMPORTANTE DO PROJETO PARA O PSO-RID
% -------------------------------------------------------------------------
%
% E O UNICO CASO QUE USA OS TRES TIPOS DE VARIAVEL AO MESMO TEMPO — o R, o I
% e o D que dao nome ao metodo:
%
%   D (diametro externo da mola) .... REAL      ->  1 dimensao real
%   N (numero de espiras) ........... INTEIRA   ->  7 bits  ([1,70])
%   d (diametro do arame) ........... DISCRETA  ->  6 bits  (42 valores)
%                                                  = 14 dimensoes
%
% [DF2011] Sec. 5.2: "the numbers of dimensional positions required for D, N
% and d are, respectively, 1 real bit, 7 binary bits, and 6 binary bits, thus
% the total number of positional dimensions of a particle becomes 14".
%
% CONSEQUENCIA PRATICA: e o PRIMEIRO caso do projeto com variavel REAL, e
% portanto o primeiro em que tres mecanismos ate agora INERTES entram em
% operacao — os mesmos que foram avaliados e mantidos em 2026-09-01 sem que
% nenhum experimento os exercitasse:
%
%   1. LIMITE DE VELOCIDADE V_max  (Shi & Eberhart 1998)
%   2. CLAMP DE POSICAO REAL       (obrigatorio: sem ele mutacao_polinomial
%                                   produz numero COMPLEXO silenciosamente)
%   3. MUTACAO POLINOMIAL Eq. (9)  com eta sorteado em [25,45]
%
% Ate aqui, todos os casos estruturais usavam so 'D' e o trem de engrenagens
% so 'I'. Este e o primeiro em que o ramo real de pso_rid.m roda de verdade.
%
% -------------------------------------------------------------------------
% O PROBLEMA ([DF2011] Sec. 5.2, Fig. 2, Eq. 13, Tabelas 5 e 6)
% -------------------------------------------------------------------------
%
% Minimizar o VOLUME DE ARAME de uma mola helicoidal de compressao que
% suporte uma carga dada sem falhar.
%
%   Determinar x = (D, N, d)
%
%   minimizar  f(x) = (pi^2/4) * D * d^2 * (N + 2)          [in^3]
%
%   sujeito a  g1 = 8*c*K*Fmax/(pi*d^2) - S            <= 0   (tensao cisalhante)
%              g2 = l - lmax                           <= 0   (comprimento livre)
%              g3 = dmin - d                           <= 0   (arame minimo)
%              g4 = (D + d) - Dmax                     <= 0   (diametro externo)
%              g5 = 3.0 - c                            <= 0   (indice de mola)
%              g6 = delta_p - delta_pm                 <= 0   (deflexao na pre-carga)
%              g7 = delta_p + (Fmax-Fp)/k
%                        + 1.05*(N+2)*d - l            <= 0
%              g8 = delta_w - (Fmax-Fp)/k              <= 0
%
%   com        c = D/d                          (indice de mola)
%              K = (4c-1)/(4c-4) + 0.615/c      (fator de Wahl)
%              k = G*d/(8*N*c^3)                (rigidez da mola)
%              delta_p = Fp/k
%              l = Fmax/k + 1.05*(N+2)*d
%
% g7 E IDENTICAMENTE NULO. Substituindo l e delta_p nas definicoes acima:
%
%   g7 = Fp/k + (Fmax-Fp)/k + 1.05(N+2)d - [Fmax/k + 1.05(N+2)d] = 0
%
% para qualquer x. E por isso que a coluna -g7 da Tabela 6 vale 0.0 em TODOS
% os sete metodos comparados, inclusive nos que sao claramente inferiores.
% A restricao e redundante; fica implementada por fidelidade a Eq. (13) e
% para que o relatorio possa reproduzir a Tabela 6 linha a linha.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%
%   [DF2011]  Datta, D.; Figueira, J.R. Applied Soft Computing 11 (2011)
%             3625-3633. -> Sec. 5.2, Fig. 2, Eq. (13), Tabelas 5 e 6
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%   [SE1998]  Shi, Y.; Eberhart, R. "A modified particle swarm optimizer",
%             IEEE ICEC 1998, p. 69-73. (origem do V_max)
%
% UNIDADES: o problema e IMPERIAL de origem (lb, in, psi) e assim fica —
% converter para SI so introduziria arredondamento sem beneficio, e a
% Tabela 6 e publicada nessas unidades.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_datta_mola                 % semente 42, 30 execucoes (protocolo DF)
%   main_datta_mola(7, 100)
%   caso = main_datta_mola('caso');
%   r = main_datta_mola;
%
% SAIDA
%   resultado : struct com .melhor_x .melhor_f .g .viavel .taxa_sucesso
%               .f_por_run .x_por_run .historicos .caso .seed

% -------------------------------------------------------------------------
% ACESSOR DO CASO
% -------------------------------------------------------------------------

if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    resultado = caso_mola();
    return;
end

if nargin < 1 || isempty(seed),   seed   = 42; end
if nargin < 2 || isempty(n_runs), n_runs = 30; end   % [DF2011]: 30 execucoes

garantir_caminhos();
rng(seed);

caso = caso_mola();

funcao_objetivo = @(x) avaliar_projeto(x, caso);

% -------------------------------------------------------------------------
% CONFIGURACAO DO PSO-RID — protocolo de [DF2011] Sec. 5
% -------------------------------------------------------------------------
% eta em [25,45] e enxame em [50,100] por execucao, max 1000 geracoes.
% Diferente do trem de engrenagens, aqui eta TEM efeito: a mutacao polinomial
% da Eq. (9) atua sobre a dimensao real D.

pso_params = struct();
pso_params.max_iter               = 1000;     % [DF2011] Sec. 5
pso_params.w                      = 1.0;      % [DF2011] Sec. 5
pso_params.c1                     = 1.0;      % [DF2011] Sec. 5
pso_params.c2                     = 2.0;      % [DF2011] Sec. 5
pso_params.pm                     = 0.15;     % [DF2011] Sec. 4.3/5
pso_params.auto_adaptativo        = true;     % [DF2011] Sec. 5
pso_params.eta_min                = 25;       % [DF2011] Sec. 5
pso_params.eta_max                = 45;       % [DF2011] Sec. 5
pso_params.decodificacao_discreta = 'proporcional';   % ver [D1] em pso_rid.m
pso_params.tol_estagnacao         = 150;
pso_params.verbose                = false;
pso_params.print_interval         = 200;

imprimir_cabecalho(caso, seed, n_runs, pso_params);

% -------------------------------------------------------------------------
% EXECUCAO MULTI-START
% -------------------------------------------------------------------------

historicos = cell(n_runs, 1);
f_por_run  = nan(n_runs, 1);
x_por_run  = nan(n_runs, 3);
enxames    = nan(n_runs, 1);
avaliacoes = nan(n_runs, 1);
melhor_f   = inf;
melhor_x   = [];

for r = 1:n_runs
    pso_params.n_particulas = randi([50, 100]);   % [DF2011] Sec. 5
    enxames(r) = pso_params.n_particulas;

    [x_r, f_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}  = hist_r;
    x_por_run(r,:) = x_r;
    avaliacoes(r)  = det_r.n_avaliacoes;

    % So conta como resultado quem terminou VIAVEL ([DEB2000] criterio 1)
    if det_r.gbest_viol <= 0
        f_por_run(r) = f_r;
        if f_r < melhor_f
            melhor_f = f_r;
            melhor_x = x_r;
        end
    end

    fprintf('  Run %2d/%d | enxame %3d | D = %.6f  N = %2d  d = %.4f | f = %.6f | viol = %.2e | aval %6d\n', ...
            r, n_runs, enxames(r), x_r(1), x_r(2), x_r(3), f_r, det_r.gbest_viol, det_r.n_avaliacoes);
end

if isempty(melhor_x)
    error('main_datta_mola:semSolucaoViavel', ...
        'Nenhuma das %d execucoes terminou VIAVEL.', n_runs);
end

% -------------------------------------------------------------------------
% RELATORIO
% -------------------------------------------------------------------------
%
% [DF2011] Sec. 5.2: "fixed values of N = 9 and d = 0.283 in. are obtained in
% each of the 30 runs of this problem. The only variation is obtained in the
% values of D, which vary in the range of [1.223041,1.223060] in."
%
% Logo o criterio de sucesso natural aqui e acertar a PARTE DISCRETA/INTEIRA
% (N = 9 e d = 0.283), que e o que o artigo diz obter em 100% das execucoes;
% o D real fica com uma faixa de tolerancia.

acertou_NI   = (x_por_run(:,2) == caso.ref_N) & (abs(x_por_run(:,3) - caso.ref_d) < 1e-9);
taxa_sucesso = 100 * nnz(acertou_NI & isfinite(f_por_run)) / n_runs;

imprimir_relatorio(melhor_x, melhor_f, caso);
imprimir_estatisticas(f_por_run, x_por_run, enxames, avaliacoes, taxa_sucesso, caso);

fig = plot_convergencia(historicos, caso.ref_f, ...
                        'PSO-RID — Mola de Compressao [DF2011] Sec. 5.2');
salvar_figura(fig, 'convergencia_datta_mola');

% -------------------------------------------------------------------------
% SAIDA ESTRUTURADA
% -------------------------------------------------------------------------

[~, viol_final, g_final] = avaliar_projeto(melhor_x, caso);

resultado.caso         = caso;
resultado.melhor_x     = melhor_x;      % [D N d]
resultado.melhor_f     = melhor_f;
resultado.g            = g_final;       % as 8 restricoes, forma bruta
resultado.violacao     = viol_final;
resultado.viavel       = viol_final <= 0;
resultado.taxa_sucesso = taxa_sucesso;
resultado.f_por_run    = f_por_run;
resultado.x_por_run    = x_por_run;
resultado.enxames      = enxames;
resultado.avaliacoes   = avaliacoes;
resultado.historicos   = historicos;
resultado.seed         = seed;

end


% #########################################################################
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################


% -------------------------------------------------------------------------
% >>> LOCAL caso_mola — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_datta_mola('caso')
% -------------------------------------------------------------------------

function caso = caso_mola()
% CASO_MOLA  Parametros da mola de compressao de [DF2011] Sec. 5.2.

caso.nome = 'Datta & Figueira (2011) - Mola Helicoidal de Compressao (Sec. 5.2)';

% -------------------------------------------------------------------------
% DADOS NUMERICOS ([DF2011] Sec. 5.2, unidades imperiais)
% -------------------------------------------------------------------------

caso.Fmax     = 1000.0;      % [lb]  carga maxima de trabalho
caso.S        = 189000.0;    % [psi] tensao cisalhante admissivel
caso.lmax     = 14.0;        % [in]  comprimento livre maximo
caso.dmin     = 0.2;         % [in]  diametro minimo do arame
caso.Dmax     = 3.0;         % [in]  diametro externo maximo da mola
caso.Fp       = 300.0;       % [lb]  forca de pre-carga
caso.delta_pm = 6.0;         % [in]  deflexao maxima sob pre-carga
caso.delta_w  = 1.25;        % [in]  deflexao da pre-carga ate a carga maxima
caso.G        = 11.5e6;      % [psi] modulo de cisalhamento

% -------------------------------------------------------------------------
% LIMITES DAS VARIAVEIS
% -------------------------------------------------------------------------
%
% [DF2011] Sec. 5.2: "Combining constraints g3(x)-g5(x), the limits of D and
% N are fixed as [0.6,3.0] in. and [1,70], respectively. Since d is a
% discrete variable with 42 allowable values given in Table 5, its integer
% limits are set automatically as [1,42]."
%
% Note que g3, g4 e g5 CONTINUAM sendo restricoes explicitas do problema: os
% limites acima sao o envelope que elas implicam, nao um substituto delas.

caso.D_min = 0.6;    caso.D_max_var = 3.0;    % [in]
caso.N_min = 1;      caso.N_max     = 70;

% -------------------------------------------------------------------------
% CATALOGO DE DIAMETROS DE ARAME ([DF2011] Tabela 5, 42 valores em polegadas)
% -------------------------------------------------------------------------

caso.catalogo_d = [ ...
    0.0090, 0.0095, 0.0104, 0.0118, 0.0128, 0.0132, 0.0140, ...
    0.0150, 0.0162, 0.0173, 0.0180, 0.0200, 0.0230, 0.0250, ...
    0.0280, 0.0320, 0.0350, 0.0410, 0.0470, 0.0540, 0.0630, ...
    0.0720, 0.0800, 0.0920, 0.1050, 0.1200, 0.1350, 0.1480, ...
    0.1620, 0.1770, 0.1920, 0.2070, 0.2250, 0.2440, 0.2630, ...
    0.2830, 0.3070, 0.3310, 0.3620, 0.3940, 0.4375, 0.5000];

assert(numel(caso.catalogo_d) == 42, 'caso_mola:catalogo', ...
    '[DF2011] Tabela 5 lista 42 diametros de arame.');

% Observacao util: g3 (d >= dmin = 0.2 in) elimina de saida os 31 primeiros
% valores da tabela. So os 11 ultimos (0.2070 a 0.5000) podem ser viaveis.
% Mesmo assim o catalogo entra INTEIRO, como em [DF2011] — filtra-lo seria
% resolver um problema diferente, e a busca precisa aprender a evitar a
% regiao inviavel por conta propria (e o que g3 mede).

% -------------------------------------------------------------------------
% VARIAVEIS DE PROJETO — UMA DE CADA TIPO (R, I, D)
% -------------------------------------------------------------------------
% Ordem conforme a Eq. (13): x = (D, N, d)

caso.config_vars(1).tipo   = 'R';    % D — diametro externo da mola
caso.config_vars(1).min    = caso.D_min;
caso.config_vars(1).max    = caso.D_max_var;

caso.config_vars(2).tipo   = 'I';    % N — numero de espiras
caso.config_vars(2).min    = caso.N_min;
caso.config_vars(2).max    = caso.N_max;

caso.config_vars(3).tipo   = 'D';    % d — diametro do arame
caso.config_vars(3).opcoes = caso.catalogo_d;

caso.nomes_vars = {'D [in] (real)', 'N (inteira)', 'd [in] (discreta)'};

% -------------------------------------------------------------------------
% SOLUCOES DE REFERENCIA ([DF2011] Tabela 6)
% -------------------------------------------------------------------------
%
%  Fonte                            D [in]     N   d [in]    f(x) [in^3]
%  ------------------------------------------------------------------------
%  Sandgren [24]                    1.180701  10   0.283     2.7995
%  Chen e Tsao [3]                  1.2287     9   0.283     2.6709
%  Wu e Chow [27]                   1.227411   9   0.283     2.6681
%  Guo et al. [12]                  1.223      9   0.283     2.659
%  Lampinen e Zelinka [17]          1.223041   9   0.283     2.65856
%  PSO de Kennedy e Eberhart [16]   1.223047   9   0.283     2.658573
%  PSO proposto ([DF2011])          1.223041   9   0.283     2.658559
%
% [DF2011] Sec. 5.2 sobre a propria linha: N = 9 e d = 0.283 em TODAS as 30
% execucoes; so D varia, na faixa [1.223041, 1.223060], com f resultante em
% [2.658559, 2.658599].
%
% Custo computacional reportado: avaliacoes de funcao entre 4784 e 98992,
% media 51696.27 e mediana 40726 nas 30 execucoes.

caso.ref_fontes = {'Sandgren [24]', 'Chen e Tsao [3]', 'Wu e Chow [27]', ...
                   'Guo et al. [12]', 'Lampinen e Zelinka [17]', ...
                   'PSO Kennedy e Eberhart [16]', 'PSO proposto [DF2011]'};
caso.ref_x = [1.180701 10 0.283; 1.2287    9 0.283; 1.227411 9 0.283; ...
              1.223    9  0.283; 1.223041  9 0.283; 1.223047 9 0.283; ...
              1.223041 9  0.283];
caso.ref_f_pub = [2.7995 2.6709 2.6681 2.659 2.65856 2.658573 2.658559];

% Melhor solucao publicada, usada como alvo
caso.ref_D = 1.223041;
caso.ref_N = 9;
caso.ref_d = 0.283;
caso.ref_f = 2.658559;
caso.ref_faixa_D = [1.223041, 1.223060];
caso.ref_faixa_f = [2.658559, 2.658599];

% -------------------------------------------------------------------------
% ESTRUTURA DO ESPACO DE BUSCA (enumeracao de 2026-09-02)
% -------------------------------------------------------------------------
%
% Como f cresce monotonicamente com D, para cada par (N,d) o melhor D e o
% MENOR viavel. Isso permite enumerar o problema inteiro: 70 valores de N x
% 42 de d, com D otimizado por bissecao. Resultado:
%
%   pares (N,d) viaveis ....................... 192 de 2940
%   otimo global .............................. N=9,  d=0.2830   f = 2.6586
%   2o lugar .................................. N=5,  d=0.3070   f = 2.7001  (+1.52%)
%   3o lugar .................................. N=10, d=0.2830   f = 2.8004  (+5.29%)
%   4o lugar .................................. N=6,  d=0.3070   f = 2.9038  (+9.18%)
%   5o lugar .................................. N=17, d=0.2630   f = 2.9109  (+9.44%)
%   pares dentro de 10% do otimo .............. 5
%
% A enumeracao CONFIRMA que N=9 e d=0.283 e o otimo global, como [DF2011]
% afirma. Mas mostra tambem por que o problema e traicoeiro: o segundo melhor
% par esta a apenas 1.52% do primeiro. Distinguir os dois exige resolver bem
% a variavel REAL D dentro de cada bacia antes de comparar as bacias.
%
% -------------------------------------------------------------------------
% DESEMPENHO OBTIDO vs [DF2011]
% -------------------------------------------------------------------------
%
% Diferente do trem de engrenagens (Sec. 5.1), aqui o solver ATINGE o otimo
% publicado: f = 2.658559 in^3, igual a Tabela 6 nas seis casas, e TODAS as
% execucoes terminam viaveis. O que nao se reproduz e a CONSISTENCIA.
%
% [DF2011] diz obter N=9 e d=0.283 em todas as 30 execucoes. Medido aqui em
% 100 execucoes (2026-09-02):
%
%   acerta N=9 e d=0.283 ......  4/100 = 4.0%   IC95% [1.6%, 9.8%]
%   termina viavel ...........  100/100 = 100%
%   f melhor .................  2.658559  (identico ao publicado)
%   f mediana ................  2.903580  (+9.2%)
%
% As demais execucoes param em bacias vizinhas — 11 pares (N,d) distintos
% apareceram, com a seguinte distribuicao:
%
%   N= 6  d=0.3070   28/100    (4o melhor par,  +9.18%)
%   N= 4  d=0.3310   15/100    (13o melhor par, +20.48%)
%   N= 5  d=0.3070   13/100    (2o melhor par,  +1.52%)
%   N=10  d=0.2830   11/100    (3o melhor par,  +5.29%)
%   N= 8  d=0.3070   10/100
%   N=17  d=0.2630    9/100    (5o melhor par,  +9.44%)
%
% ATENCAO A ESTE PONTO: o par mais frequente nao e o segundo melhor, e o
% QUARTO; e 15% das execucoes param num par 20% pior. Ou seja, o problema nao
% e apenas "quase-empate dificil de desempatar" — o enxame nem sequer se
% concentra nas melhores bacias. Uma amostra menor engana: com 30 execucoes a
% taxa medida foi 13.3%, tres vezes maior que a de 100 execucoes, e dentro do
% IC95% de ambas. Ao comparar configuracoes neste caso, use amostras grandes.
%
% O custo computacional, esse sim, e comparavel ao do artigo:
%
%                          avaliacoes: faixa        media      mediana
%   [DF2011] Sec. 5.2      [4784, 98992]            51696.27   40726
%   este solver (30 runs)  [9424, 64228]            40152.33   43626
%
% Ou seja, gasta-se o mesmo esforco e chega-se ao mesmo melhor valor, mas com
% muito mais dispersao entre execucoes. E o mesmo padrao ja registrado em
% main_datta_engrenagens.m, so que bem menos severo — la o otimo aparecia em
% 1% das execucoes, aqui em 13%. Por isso o relatorio deste orquestrador
% imprime a TAXA DE SUCESSO e a distribuicao de pares (N,d): o numero
% interessante nao e o melhor f, e com que frequencia ele reaparece.

end


% -------------------------------------------------------------------------
% >>> LOCAL restricoes — Eq. (13) de [DF2011], forma BRUTA
% -------------------------------------------------------------------------

function [g, aux] = restricoes(x, caso)
% RESTRICOES  As oito g(x) da Eq. (13), nas unidades originais do artigo.
%
% Devolve tambem as grandezas derivadas, que o relatorio usa para reproduzir
% a Tabela 6.

D = x(1);  N = x(2);  d = x(3);

aux.c       = D / d;                                    % indice de mola
aux.K       = (4*aux.c - 1)/(4*aux.c - 4) + 0.615/aux.c; % fator de Wahl
aux.k       = caso.G * d / (8 * N * aux.c^3);            % rigidez [lb/in]
aux.delta_p = caso.Fp / aux.k;                           % [in]
aux.l       = caso.Fmax / aux.k + 1.05*(N + 2)*d;        % comprimento livre [in]
aux.tau     = 8 * aux.c * aux.K * caso.Fmax / (pi * d^2);% tensao cisalhante [psi]

g = zeros(1, 8);
g(1) = aux.tau - caso.S;                                       % [psi]
g(2) = aux.l - caso.lmax;                                      % [in]
g(3) = caso.dmin - d;                                          % [in]
g(4) = (D + d) - caso.Dmax;                                    % [in]
g(5) = 3.0 - aux.c;                                            % [-]
g(6) = aux.delta_p - caso.delta_pm;                            % [in]
g(7) = aux.delta_p + (caso.Fmax - caso.Fp)/aux.k ...
       + 1.05*(N + 2)*d - aux.l;                               % [in] == 0
g(8) = caso.delta_w - (caso.Fmax - caso.Fp)/aux.k;             % [in]
end


% -------------------------------------------------------------------------
% >>> LOCAL avaliar_projeto — funcao objetivo
% -------------------------------------------------------------------------

function [custo, violacao, g] = avaliar_projeto(x, caso)
% AVALIAR_PROJETO  Volume de arame e violacao das oito restricoes.
%
%   f(x) = (pi^2/4) * D * d^2 * (N + 2)      [in^3]     [DF2011] Eq. (13)
%
% -------------------------------------------------------------------------
% MEDIDA DE VIOLACAO: NORMALIZADA, e nao a soma bruta das g(x)
% -------------------------------------------------------------------------
% As oito restricoes da Eq. (13) estao em unidades DIFERENTES: g1 em psi
% (ordem de 1e5), g2/g3/g4/g6/g7/g8 em polegadas (ordem de 1e0) e g5
% adimensional. Somar os excessos crus faria g1 dominar por cinco ordens de
% grandeza: 1 psi de excesso de tensao pesaria o mesmo que 1 polegada de
% excesso de comprimento, embora a primeira consuma 0.0005% da folga e a
% segunda 7%.
%
% Isso importa porque o criterio 3 de [DEB2000] ordena solucoes inviaveis
% SOMENTE pela violacao. E exatamente o defeito corrigido em [D8] para os
% casos estruturais (ver avaliar_projeto de main_hadi_nao_linear.m); aqui ele
% seria muito pior, porque la havia duas familias de unidades e aqui ha tres.
%
% Cada g e entao reescrita na forma adimensional (grandeza/limite - 1):
%
%   g1 -> tau/S - 1                     g5 -> 3.0/c - 1
%   g2 -> l/lmax - 1                    g6 -> delta_p/delta_pm - 1
%   g3 -> dmin/d - 1                    g7 -> identicamente 0 (ver cabecalho)
%   g4 -> (D+d)/Dmax - 1                g8 -> delta_w*k/(Fmax-Fp) - 1
%
% A viabilidade (o SINAL de cada g) e identica nas duas formas — todos os
% denominadores sao constantes positivas do problema. So a MAGNITUDE muda, e
% e ela que o criterio 3 usa para ordenar.

D = x(1);  N = x(2);  d = x(3);

custo = (pi^2 / 4) * D * d^2 * (N + 2);

[g, aux] = restricoes(x, caso);

gn = zeros(1, 8);
gn(1) = aux.tau / caso.S - 1;
gn(2) = aux.l   / caso.lmax - 1;
gn(3) = caso.dmin / d - 1;
gn(4) = (D + d) / caso.Dmax - 1;
gn(5) = 3.0 / aux.c - 1;
gn(6) = aux.delta_p / caso.delta_pm - 1;
gn(7) = g(7) / caso.lmax;                 % nulo por construcao; normalizado
gn(8) = caso.delta_w * aux.k / (caso.Fmax - caso.Fp) - 1;

violacao = sum(max(0, gn));
end


% -------------------------------------------------------------------------
% >>> LOCAIS de saida
% -------------------------------------------------------------------------

function imprimir_cabecalho(caso, seed, n_runs, pso_params)
aux = pso_rid('auxiliares');
[mapa, total_dim] = aux.mapear_dimensoes(caso.config_vars);
fprintf('\n============================================================\n');
fprintf(' %s\n', caso.nome);
fprintf('============================================================\n');
fprintf(' Objetivo           : min (pi^2/4)*D*d^2*(N+2)  [in^3]   (Eq. 13)\n');
fprintf(' Restricoes         : 8 (g7 e identicamente nula — ver cabecalho)\n');
fprintf(' Variaveis          : os TRES tipos do PSO-RID\n');
for k = 1:numel(caso.config_vars)
    if strcmp(mapa(k).tipo, 'R')
        fprintf('    %-18s REAL      [%g, %g]        -> 1 dim real\n', ...
                caso.nomes_vars{k}, caso.config_vars(k).min, caso.config_vars(k).max);
    elseif strcmp(mapa(k).tipo, 'I')
        fprintf('    %-18s INTEIRA   [%d, %d]         -> %d bits\n', ...
                caso.nomes_vars{k}, caso.config_vars(k).min, ...
                caso.config_vars(k).max, mapa(k).n_bits);
    else
        fprintf('    %-18s DISCRETA  %d opcoes        -> %d bits\n', ...
                caso.nomes_vars{k}, numel(caso.config_vars(k).opcoes), mapa(k).n_bits);
    end
end
fprintf(' Total              : %d dimensoes   ([DF2011] Sec. 5.2 diz 14)\n', total_dim);
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d   (protocolo de [DF2011] Sec. 5)\n', n_runs);
fprintf(' Enxame             : sorteado em [50,100] a cada execucao\n');
fprintf(' Iteracoes maximas  : %d | estagnacao: %d\n', ...
        pso_params.max_iter, pso_params.tol_estagnacao);
fprintf(' Melhor publicado   : D=%.6f N=%d d=%.4f -> f=%.6f in^3\n', ...
        caso.ref_D, caso.ref_N, caso.ref_d, caso.ref_f);
fprintf('============================================================\n\n');
end


function imprimir_relatorio(x, f, caso)
fprintf('\n============================================================\n');
fprintf(' RESULTADO — comparacao com a Tabela 6 de [DF2011]\n');
fprintf('============================================================\n');
fprintf(' %-30s %10s %4s %8s %12s\n', 'Fonte', 'D [in]', 'N', 'd [in]', 'f [in^3]');
fprintf(' --------------------------------------------------------------------\n');
for i = 1:numel(caso.ref_fontes)
    xi = caso.ref_x(i,:);
    fi = (pi^2/4) * xi(1) * xi(3)^2 * (xi(2) + 2);
    fprintf(' %-30s %10.6f %4d %8.4f %12.6f\n', caso.ref_fontes{i}, xi(1), xi(2), xi(3), fi);
end
fprintf(' --------------------------------------------------------------------\n');
fprintf(' %-30s %10.6f %4d %8.4f %12.6f\n', '>> ESTE PSO-RID', x(1), x(2), x(3), f);
fprintf(' --------------------------------------------------------------------\n');
fprintf(' Faixa publicada em 30 execucoes: D em [%.6f, %.6f], f em [%.6f, %.6f]\n', ...
        caso.ref_faixa_D, caso.ref_faixa_f);
fprintf(' Diferenca do melhor publicado  : %+.6f in^3 (%+.4f%%)\n', ...
        f - caso.ref_f, 100*(f - caso.ref_f)/caso.ref_f);

% Reproducao da Tabela 6: valores de -g(x)
[g, ~] = restricoes(x, caso);
fprintf('\n RESTRICOES  -g(x)  (comparar com a Tabela 6; >= 0 = satisfeita)\n');
fprintf(' --------------------------------------------------------------------\n');
rot = {'g1 tensao cisalhante [psi]', 'g2 comprimento livre  [in]', ...
       'g3 arame minimo       [in]', 'g4 diametro externo   [in]', ...
       'g5 indice de mola      [-]', 'g6 defl. pre-carga    [in]', ...
       'g7 (redundante)       [in]', 'g8 defl. de trabalho  [in]'};
for i = 1:8
    fprintf('   %-28s %14.6f  %s\n', rot{i}, -g(i), ...
            ternario(g(i) <= 1e-9, 'ok', '<< VIOLADA'));
end
fprintf('============================================================\n\n');
end


function s = ternario(cond, a, b)
if cond, s = a; else, s = b; end
end


function imprimir_estatisticas(f_por_run, x_por_run, enxames, avaliacoes, taxa, caso)
n = numel(f_por_run);
val = f_por_run(isfinite(f_por_run));
fprintf(' ESTATISTICAS DAS %d EXECUCOES\n', n);
fprintf(' ------------------------------------------------------------\n');
fprintf('  Execucoes viaveis : %d de %d\n', numel(val), n);
fprintf('  TAXA DE SUCESSO   : %.1f%% acertaram N=%d e d=%.4f\n', taxa, caso.ref_N, caso.ref_d);
fprintf('    [DF2011] Sec. 5.2 reporta 100%% (N e d fixos em 30/30)\n');
fprintf('  f melhor          : %.6f in^3   (publicado %.6f)\n', min(val), caso.ref_f);
fprintf('  f mediana         : %.6f\n', median(val));
fprintf('  f media           : %.6f\n', mean(val));
fprintf('  f desvio padrao   : %.3e\n', std(val));
fprintf('  f pior            : %.6f\n', max(val));
fprintf('  D                 : [%.6f, %.6f]   (publicado [%.6f, %.6f])\n', ...
        min(x_por_run(:,1)), max(x_por_run(:,1)), caso.ref_faixa_D);
fprintf('  Enxame            : %d a %d particulas\n', min(enxames), max(enxames));
fprintf('  Avaliacoes        : %d a %d | media %.2f | mediana %.0f\n', ...
        min(avaliacoes), max(avaliacoes), mean(avaliacoes), median(avaliacoes));
fprintf('    [DF2011] reporta [4784, 98992], media 51696.27, mediana 40726\n');

[u, ~, idx] = unique(x_por_run(:,2:3), 'rows');
fprintf('  Pares (N,d) distintos encontrados: %d\n', size(u,1));
for i = 1:size(u,1)
    fprintf('     N=%2d  d=%.4f  em %2d/%d execucoes\n', u(i,1), u(i,2), nnz(idx==i), n);
end
fprintf(' ------------------------------------------------------------\n\n');
end


function garantir_caminhos()
% Este orquestrador precisa apenas do Bloco 1 (nao usa o FEM).
if exist('pso_rid', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
