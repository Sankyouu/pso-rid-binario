function resultado = main_datta_engrenagens(seed, n_runs)
% MAIN_DATTA_ENGRENAGENS  Projeto do trem de engrenagens de [DF2011] Sec. 5.1.
%
% -------------------------------------------------------------------------
% O QUE ESTE ARQUIVO FAZ
% -------------------------------------------------------------------------
%
% A funcao objetivo e analitica — uma razao de engrenagens — entao nao ha solver
% estrutural envolvido, e todo o experimento cabe neste unico arquivo:
%
%   Bloco 1 (otimizador) : pso_rid.m   <- 01_pso_rid/
%   Funcao objetivo      : local, aqui embaixo (avaliar_projeto)
%   Definicao do problema: local, aqui embaixo (caso_engrenagens)
%
% -------------------------------------------------------------------------
% POR QUE ESTE CASO IMPORTA PARA O PROJETO
% -------------------------------------------------------------------------
%
% 1. E O PRIMEIRO CASO COM VARIAVEL INTEIRA ('I'). Todos os outros
%    (Hadi 10, Hadi 20, Hadi 51, Awruch) usam so variaveis discretas 'D'.
%    O tipo 'I' existia no solver, era exercitado apenas por teste unitario,
%    e ficava inerte em todos os experimentos reais. Aqui ele e o problema.
%
% 2. E O PRIMEIRO CASO EM QUE viol_estrutural REALMENTE DISPARA. Com
%    12 <= z <= 60 e 6 bits por variavel, a Eq. (6) de [DF2011] gera valores
%    de 0 a 63: os codigos 0..11 e 61..63 caem FORA do dominio. Pela regra
%    [D2], esses valores nao sao saturados — sao sinalizados como inviaveis e
%    a funcao objetivo nem chega a ser chamada, coerente com [DEB2000] pag.
%    316 e com a Eq. (12) abaixo, que trata os limites como RESTRICAO g(x).
%
%    Fracao de codigos validos por variavel : 49/64 = 76.6%
%    Fracao de particulas inteiramente valida: 0.766^4 = 34.4%
%
%    Ou seja, em torno de dois tercos das particulas caem em custo Inf. Isso
%    e consequencia direta da formulacao de [DF2011], nao uma escolha nossa.
%
% 3. E O CASO COM VERIFICACAO MAIS FORTE DO PROJETO. O otimo global e
%    conhecido e provado ([DF2011] cita [10]); o artigo reporta atingi-lo em
%    TODAS as 30 execucoes. Logo a metrica aqui nao e "que peso deu", e sim
%    TAXA DE SUCESSO em 30 execucoes — comparavel de forma direta.
%
% -------------------------------------------------------------------------
% O PROBLEMA ([DF2011] Sec. 5.1, Fig. 1, Eqs. 11-12)
% -------------------------------------------------------------------------
%
% Trem de engrenagens composto, dois pares de rodas: a-c e b-d, com a e b
% MOTORAS e c e d MOVIDAS.
%
%   Eq. (11):  razao = (dentes das motoras)/(dentes das movidas) = za*zb/(zc*zd)
%
%   Eq. (12):  determinar x = (za, zb, zc, zd)
%              minimizar  f(x) = [ 1/6.931 - za*zb/(zc*zd) ]^2
%              sujeito a  g(x) == 12 <= za, zb, zc, zd <= 60
%                         za, zb, zc, zd inteiros
%
% "It is required to determine the values of za, zb, zc and zd so that a gear
%  ratio, as close as possible to 1/6.931, can be obtained. The only
%  constraint in the problem is that the number of teeth on any gear should
%  be in the range of [12,60]."
%
% Note que f e um ERRO QUADRATICO, nao uma restricao: nao ha violacao a
% reportar. A unica restricao do problema sao os limites inteiros, e eles ja
% sao tratados dentro do decodificador. Por isso avaliar_projeto devolve
% violacao identicamente nula — ver a nota la embaixo.
%
% -------------------------------------------------------------------------
% REFERENCIAS
% -------------------------------------------------------------------------
%
%   [DF2011]  Datta, D.; Figueira, J.R. "A real-integer-discrete-coded
%             particle swarm optimization for design problems". Applied Soft
%             Computing 11 (2011) 3625-3633.  ->  Sec. 5.1, Fig. 1, Tabela 4
%             -> 00_docs/artigos/Areal-integer-discrete-coded particle swarm
%                optimization for design problems.pdf
%   [DEB2000] Deb, K., Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_datta_engrenagens              % semente 42, 30 execucoes (protocolo
%                                       % do artigo)
%   main_datta_engrenagens(7, 100)      % outra semente, 100 execucoes
%   caso = main_datta_engrenagens('caso');
%   r = main_datta_engrenagens;
%
% SAIDA
%   resultado : struct com .melhor_z .melhor_f .razao .erro_pct .taxa_sucesso
%               .f_por_run .z_por_run .historicos .caso .seed

% -------------------------------------------------------------------------
% ACESSOR DO CASO
% -------------------------------------------------------------------------

if nargin == 1 && ischar(seed) && strcmp(seed, 'caso')
    resultado = caso_engrenagens();
    return;
end

if nargin < 1 || isempty(seed),   seed   = 42; end
if nargin < 2 || isempty(n_runs), n_runs = 30; end   % [DF2011]: 30 execucoes

garantir_caminhos();
rng(seed);

caso = caso_engrenagens();

funcao_objetivo = @(z) avaliar_projeto(z, caso);

% -------------------------------------------------------------------------
% CONFIGURACAO DO PSO-RID — protocolo de [DF2011] Sec. 5
% -------------------------------------------------------------------------
%
% "The distribution index in Eq. (9) is assigned a random value, in the range
%  of [25,45], in different runs of a problem. Similarly, the PSO swarm size
%  in different runs is randomly fixed in the range of [50,100]. In the case
%  of execution time, each run is allowed to be continued for a maximum
%  number of 1000 generations. However, a run is terminated in between when
%  no improvement in the best objective value is noticed."
%
% O tamanho do enxame e sorteado EM CADA EXECUCAO, dentro do laco abaixo.
% O indice de distribuicao eta ja e sorteado internamente por pso_rid, nos
% limites eta_min/eta_max (padrao 25/45, que sao exatamente os do artigo) —
% aqui ele nao tem efeito pratico, porque a mutacao polinomial da Eq. (9) so
% se aplica a dimensoes REAIS e este problema nao tem nenhuma.

pso_params = struct();
pso_params.max_iter               = 1000;     % [DF2011] Sec. 5
pso_params.w                      = 1.0;      % [DF2011] Sec. 5
pso_params.c1                     = 1.0;      % [DF2011] Sec. 5
pso_params.c2                     = 2.0;      % [DF2011] Sec. 5
pso_params.pm                     = 0.15;     % [DF2011] Sec. 4.3/5
pso_params.auto_adaptativo        = true;     % [DF2011] Sec. 5
pso_params.eta_min                = 25;       % [DF2011] Sec. 5
pso_params.eta_max                = 45;       % [DF2011] Sec. 5
pso_params.tol_estagnacao         = 100;      % ver nota em imprimir_cabecalho
pso_params.verbose                = false;    % 30 execucoes: log por execucao
pso_params.print_interval         = 200;

imprimir_cabecalho(caso, seed, n_runs, pso_params);

% -------------------------------------------------------------------------
% EXECUCAO MULTI-START
% -------------------------------------------------------------------------

historicos = cell(n_runs, 1);
f_por_run  = nan(n_runs, 1);
z_por_run  = nan(n_runs, 4);
enxames    = nan(n_runs, 1);
melhor_f   = inf;
melhor_z   = [];

for r = 1:n_runs
    % [DF2011] Sec. 5: tamanho do enxame sorteado em [50,100] a cada execucao
    pso_params.n_particulas = randi([50, 100]);
    enxames(r) = pso_params.n_particulas;

    [z_r, f_r, hist_r, det_r] = pso_rid(funcao_objetivo, caso.config_vars, pso_params);

    historicos{r}  = hist_r;
    f_por_run(r)   = f_r;
    z_por_run(r,:) = z_r;

    if f_r < melhor_f
        melhor_f = f_r;
        melhor_z = z_r;
    end

    fprintf('  Run %2d/%d | enxame %3d | z = [%2d %2d %2d %2d] | f = %.4e | razao = %.6f | iters %4d | overflow %d\n', ...
            r, n_runs, enxames(r), z_r(1), z_r(2), z_r(3), z_r(4), f_r, ...
            razao(z_r), det_r.iter_executadas, det_r.n_overflow);
end

% -------------------------------------------------------------------------
% RELATORIO
% -------------------------------------------------------------------------

% Sucesso = atingiu o otimo global conhecido (comparacao pelo VALOR de f,
% nao pelo vetor z: [DF2011] Tabela 4 mostra que (19,16,...) e (16,19,...)
% dao o mesmo f, pois za*zb e zc*zd sao simetricos em seus pares).
sucesso       = f_por_run <= caso.f_otimo * (1 + 1e-9);
taxa_sucesso  = 100 * nnz(sucesso) / n_runs;

imprimir_relatorio(melhor_z, melhor_f, caso);
imprimir_estatisticas(f_por_run, z_por_run, enxames, taxa_sucesso, caso);

fig = plot_convergencia_log(historicos, caso.f_otimo);
salvar_figura(fig, 'convergencia_datta_engrenagens');

% -------------------------------------------------------------------------
% SAIDA ESTRUTURADA
% -------------------------------------------------------------------------

resultado.caso          = caso;
resultado.melhor_z      = melhor_z;
resultado.melhor_f      = melhor_f;
resultado.razao         = razao(melhor_z);
resultado.erro_pct      = 100 * abs(razao(melhor_z) - caso.razao_alvo) / caso.razao_alvo;
resultado.taxa_sucesso  = taxa_sucesso;
resultado.f_por_run     = f_por_run;
resultado.z_por_run     = z_por_run;
resultado.enxames       = enxames;
resultado.historicos    = historicos;
resultado.seed          = seed;

end


% #########################################################################
% FUNCOES LOCAIS DO ORQUESTRADOR
% #########################################################################


% -------------------------------------------------------------------------
% >>> LOCAL caso_engrenagens — DEFINICAO DO PROBLEMA
%     Acesso externo: caso = main_datta_engrenagens('caso')
% -------------------------------------------------------------------------

function caso = caso_engrenagens()
% CASO_ENGRENAGENS  Parametros do trem de engrenagens de [DF2011] Sec. 5.1.
%
% ARQUIVO DE PARAMETROS. Sem geometria, sem material, sem malha: o problema
% e puramente numerico.
%
% SAIDA
%   caso : struct com .razao_alvo .z_min .z_max .config_vars .ref_* .f_otimo

caso.nome = 'Datta & Figueira (2011) - Trem de Engrenagens (Sec. 5.1)';

% -------------------------------------------------------------------------
% ALVO E LIMITES ([DF2011] Eq. 12)
% -------------------------------------------------------------------------

caso.razao_alvo = 1/6.931;   % razao de engrenagens desejada
caso.z_min      = 12;        % dentes minimos em qualquer roda
caso.z_max      = 60;        % dentes maximos em qualquer roda

caso.nomes_vars = {'za (motora)', 'zb (motora)', 'zc (movida)', 'zd (movida)'};

% -------------------------------------------------------------------------
% VARIAVEIS DE PROJETO — 4 INTEIRAS
% -------------------------------------------------------------------------
%
% [DF2011] Sec. 5.1: "As per the formulation of the proposed PSO, 6 binary
% bits are required to represent each variable, thus the total number of
% positional dimensions of a particle becomes 24 to represent all the four
% variables."
%
% Confere com mapear_dimensoes: para tipo 'I', n_bits = ceil(log2(max+1)) =
% ceil(log2(61)) = 6, e 4 x 6 = 24 dimensoes. Este caso e justamente o
% exemplo numerico que documenta aquela formula em pso_rid.m.

caso.config_vars = struct('tipo', {}, 'min', {}, 'max', {});
for k = 1:4
    caso.config_vars(k).tipo = 'I';
    caso.config_vars(k).min  = caso.z_min;
    caso.config_vars(k).max  = caso.z_max;
end

% -------------------------------------------------------------------------
% SOLUCOES DE REFERENCIA ([DF2011] Tabela 4)
% -------------------------------------------------------------------------
%
%  Fonte                          za  zb  zc  zd    f(x)        razao      erro
%  ----------------------------------------------------------------------------
%  Sandgren [24]                  18  22  45  60    5.7e-06     0.146666   1.65%
%  Loh e Papalambros [22]         19  16  42  50    0.23e-06    0.144762   0.334%
%  Zhang e Wang [29]              30  15  52  60    2.4e-09     0.144231   0.033%
%  Lin et al. [20]                19  16  43  49    2.7e-12     0.144281   0.0011%
%  Guo et al. [12]                16  19  43  49    2.7e-12     0.144281   0.0011%
%  PSO de Kennedy e Eberhart [16] 16  19  43  49    2.7e-12     0.144281   0.0011%
%  PSO proposto ([DF2011])        16  19  43  49    2.7e-12     0.144281   0.0011%
%
% As quatro ultimas linhas sao a MESMA solucao: (19,16,43,49) e (16,19,43,49)
% diferem so pela ordem dentro do par de motoras, e f depende apenas dos
% PRODUTOS za*zb e zc*zd.
%
% [DF2011] Sec. 5.1 sobre esta linha: "the known best solution of the problem
% could be obtained in all of its 30 runs... In fact, no further improvement
% is possible as this known best solution is the global optimum of the
% problem [10]."
%
% E ESSE O PONTO DE COMPARACAO DESTE ORQUESTRADOR: nao basta chegar perto,
% tem que acertar o otimo global, e a metrica e em quantas das 30 execucoes.

caso.ref_fontes = {'Sandgren [24]', 'Loh e Papalambros [22]', 'Zhang e Wang [29]', ...
                   'Lin et al. [20]', 'Guo et al. [12]', ...
                   'PSO Kennedy e Eberhart [16]', 'PSO proposto [DF2011]'};
caso.ref_z = [18 22 45 60; 19 16 42 50; 30 15 52 60; ...
              19 16 43 49; 16 19 43 49; 16 19 43 49; 16 19 43 49];

% Otimo global. Calculado da propria definicao em vez de copiar o "2.7e-12"
% arredondado da Tabela 4, para que a comparacao de sucesso seja exata.
caso.z_otimo = [16 19 43 49];
caso.f_otimo = (caso.razao_alvo - (16*19)/(43*49))^2;

% -------------------------------------------------------------------------
% LACUNA DE REPRODUTIBILIDADE: 1% DE SUCESSO AQUI CONTRA 100% EM [DF2011]
% (medicao de 2026-09-02; 100 execucoes por configuracao, IC de Wilson 95%)
% -------------------------------------------------------------------------
%
% [DF2011] Sec. 5.1 afirma atingir o otimo global em TODAS as 30 execucoes.
% Este solver atinge em 1 de 100:
%
%   PSO-RID (padrao do projeto) ....  1/100 =  1.0%   IC95% [0.2%,  5.4%]
%   busca aleatoria, mesmo orcamento  0/100 =  0.0%   (0.3% analitico)
%   [DF2011] reportado .............. 30/30 = 100%
%
% O PSO fica cerca de 3x acima do acaso — ele faz alguma coisa — mas nao
% chega perto do reportado. Vale entender a dificuldade do problema:
%
%   O otimo global e UNICO. Enumerando os 49^4 = 5 764 801 pontos viaveis,
%   so (16,19,43,49) atinge f = 2.70e-12, e as 4 permutacoes dentro dos pares
%   (za<->zb, zc<->zd) dao o mesmo valor. Sao 4 pontos em 5.76 milhoes.
%   Com ~4800 avaliacoes por execucao, acertar por sorte tem 0.3% de chance.
%
% O QUE FOI TESTADO E NAO EXPLICA A DIFERENCA:
%
%   - PARADA ANTECIPADA. Rodando as 1000 geracoes inteiras (sem parada por
%     estagnacao) a taxa nao melhora, e o custo sobe 6x (30664 contra 4784
%     avaliacoes por execucao).
%   - REINICIALIZACAO DE PARTICULAS ([D6]). Com e sem: 1/100 e 2/100, com
%     IC95% amplamente sobrepostos. Nao ha efeito detectavel. (Primeira
%     evidencia empirica sobre [D6], que seguia sem medicao.)
%   - AUTO-ADAPTATIVO. Desligar PIORA muito (0% em 30 execucoes), o que
%     confirma que o esquema de [DF2011] Sec. 5 ajuda.
%   - p_busca_local e pm. Varrendo 0.05/0.15/0.40 e pm 0.05/0.15: nenhuma
%     configuracao sai da faixa de poucos por cento.
%   - TRATAMENTO DO LIMITE. Reformular g(x) = 12 <= z <= 60 como restricao de
%     [DEB2000] (dominio inteiro 0..63 e violacao proporcional ao quanto
%     passou) em vez de dominio da variavel (estouro -> Inf) da 0/100. A
%     hipotese de que o Inf sufoca o gradiente rumo ao dominio viavel foi
%     TESTADA E REJEITADA para este problema.
%   - BUSCA LOCAL EM DENTES. Aplicar hill-climb de +-1 dente sobre o
%     resultado de cada execucao nao muda nada (1/100 -> 1/100): as solucoes
%     encontradas nao estao na vizinhanca do otimo, estao em outra parte da
%     variedade quase-otima.
%
% DIAGNOSTICO. A funcao objetivo tem uma variedade enorme de combinacoes
% distintas com f entre 1e-9 e 1e-10 espalhadas pelo dominio; o enxame
% converge para uma delas e nao existe caminho local ate o otimo global. Um
% detalhe agrava: no encaixe binario, mudar UM dente custa 1.96 bits em
% media e ate 6 (31->32 inverte os seis), entao a vizinhanca de bits que a
% busca local de [DF2011] explora nao corresponde a vizinhanca de dentes.
%
% NAO FOI POSSIVEL FECHAR a diferenca com nenhum ajuste dentro do que
% [DF2011] descreve. Fica registrado como divergencia em aberto. O
% orquestrador reporta TAXA DE SUCESSO justamente para deixar isso visivel a
% cada execucao, em vez de esconder atras de uma media de f.

end


% -------------------------------------------------------------------------
% >>> LOCAL avaliar_projeto — funcao objetivo, [DF2011] Eq. (12)
% -------------------------------------------------------------------------

function [custo, violacao] = avaliar_projeto(z, caso)
% AVALIAR_PROJETO  Erro quadratico da razao de engrenagens.
%
%   f(x) = [ 1/6.931 - za*zb/(zc*zd) ]^2         [DF2011] Eq. (12)
%
% VIOLACAO IDENTICAMENTE NULA, e isso NAO e um descuido. A unica restricao do
% problema e g(x) == 12 <= z <= 60 ([DF2011] Eq. 12), e ela e resolvida ANTES
% desta funcao: decodificar sinaliza qualquer dente fora de [12,60] via
% viol_estrutural, e o laco principal de pso_rid entao atribui custo = Inf e
% violacao = Inf SEM chamar a funcao objetivo — [DEB2000] pag. 316, "It does
% not make sense to compute the objective function value of an infeasible
% solution". Ver [D2] em pso_rid.m.
%
% Logo, se a execucao chegou aqui, o projeto ja e viavel por construcao.
% Devolver 0 mantem o contrato de duas saidas que domina_deb espera, sem
% inventar uma segunda medida de qualidade que o problema nao tem.

custo    = (caso.razao_alvo - (z(1)*z(2)) / (z(3)*z(4)))^2;
violacao = 0;
end


% -------------------------------------------------------------------------
% >>> LOCAL razao — Eq. (11) de [DF2011]
% -------------------------------------------------------------------------

function r = razao(z)
% Razao do trem: dentes das MOTORAS (a,b) sobre dentes das MOVIDAS (c,d).
r = (z(1)*z(2)) / (z(3)*z(4));
end


% -------------------------------------------------------------------------
% >>> LOCAIS de saida
% -------------------------------------------------------------------------

function imprimir_cabecalho(caso, seed, n_runs, pso_params)
n_bits = max(ceil(log2(caso.z_max + 1)), 1);
fprintf('\n============================================================\n');
fprintf(' %s\n', caso.nome);
fprintf('============================================================\n');
fprintf(' Objetivo           : min [1/6.931 - za*zb/(zc*zd)]^2   (Eq. 12)\n');
fprintf(' Razao alvo         : 1/6.931 = %.9f\n', caso.razao_alvo);
fprintf(' Variaveis          : 4 INTEIRAS em [%d, %d]\n', caso.z_min, caso.z_max);
fprintf(' Codificacao        : %d bits/variavel -> %d dimensoes (Eq. 5-6)\n', ...
        n_bits, 4*n_bits);
fprintf(' Codigos validos    : %d de %d por variavel (%.1f%%)\n', ...
        caso.z_max - caso.z_min + 1, 2^n_bits, ...
        100*(caso.z_max-caso.z_min+1)/2^n_bits);
fprintf(' Particula 100%% valida: %.1f%% das amostras (os demais -> custo Inf)\n', ...
        100*((caso.z_max-caso.z_min+1)/2^n_bits)^4);
fprintf(' Otimo global       : z = [%d %d %d %d], f = %.6e\n', ...
        caso.z_otimo, caso.f_otimo);
fprintf(' Semente (rng)      : %d\n', seed);
fprintf(' Execucoes          : %d   (protocolo de [DF2011] Sec. 5)\n', n_runs);
fprintf(' Enxame             : sorteado em [50,100] a cada execucao\n');
fprintf(' Iteracoes maximas  : %d | parada por estagnacao: %d iteracoes\n', ...
        pso_params.max_iter, pso_params.tol_estagnacao);
fprintf('   (o artigo diz apenas "terminated when no improvement is\n');
fprintf('    noticed", sem numero; 100 e a escolha deste orquestrador)\n');
fprintf('============================================================\n\n');
end


function imprimir_relatorio(z, f, caso)
fprintf('\n============================================================\n');
fprintf(' RESULTADO — comparacao com a Tabela 4 de [DF2011]\n');
fprintf('============================================================\n');
fprintf(' %-30s %3s %3s %3s %3s  %11s  %9s  %8s\n', ...
        'Fonte', 'za', 'zb', 'zc', 'zd', 'f(x)', 'razao', 'erro %');
fprintf(' ------------------------------------------------------------------------------\n');
for i = 1:numel(caso.ref_fontes)
    zi = caso.ref_z(i,:);
    fi = (caso.razao_alvo - razao(zi))^2;
    fprintf(' %-30s %3d %3d %3d %3d  %11.3e  %9.6f  %8.4f\n', ...
            caso.ref_fontes{i}, zi, fi, razao(zi), ...
            100*abs(razao(zi)-caso.razao_alvo)/caso.razao_alvo);
end
fprintf(' ------------------------------------------------------------------------------\n');
fprintf(' %-30s %3d %3d %3d %3d  %11.3e  %9.6f  %8.4f\n', ...
        '>> ESTE PSO-RID', z, f, razao(z), ...
        100*abs(razao(z)-caso.razao_alvo)/caso.razao_alvo);
fprintf(' ------------------------------------------------------------------------------\n');
if f <= caso.f_otimo * (1 + 1e-9)
    fprintf(' OTIMO GLOBAL ATINGIDO.\n');
else
    fprintf(' Acima do otimo global por um fator de %.3g.\n', f / caso.f_otimo);
end
fprintf('============================================================\n\n');
end


function imprimir_estatisticas(f_por_run, z_por_run, enxames, taxa, caso)
fprintf(' ESTATISTICAS DAS %d EXECUCOES\n', numel(f_por_run));
fprintf(' ------------------------------------------------------------\n');
fprintf('  TAXA DE SUCESSO   : %.1f%% (%d de %d atingiram o otimo global)\n', ...
        taxa, round(taxa*numel(f_por_run)/100), numel(f_por_run));
fprintf('    [DF2011] Sec. 5.1 reporta 100%% em 30 execucoes\n');
fprintf('  f melhor          : %.6e\n', min(f_por_run));
fprintf('  f mediana         : %.6e\n', median(f_por_run));
fprintf('  f media           : %.6e\n', mean(f_por_run));
fprintf('  f desvio padrao   : %.6e\n', std(f_por_run));
fprintf('  f pior            : %.6e\n', max(f_por_run));
fprintf('  Enxame            : %d a %d particulas (sorteado)\n', ...
        min(enxames), max(enxames));

% Solucoes distintas encontradas — [DF2011] Sec. 5 comenta que "distinct sets
% of variable values are obtained in different runs" em alguns problemas.
[unicas, ~, idx] = unique(z_por_run, 'rows');
fprintf('  Solucoes distintas: %d\n', size(unicas,1));
for i = 1:size(unicas,1)
    n = nnz(idx == i);
    fi = (caso.razao_alvo - razao(unicas(i,:)))^2;
    fprintf('     [%2d %2d %2d %2d]  f = %.3e  em %2d/%d execucoes\n', ...
            unicas(i,:), fi, n, numel(f_por_run));
end
fprintf(' ------------------------------------------------------------\n\n');
end


function fig = plot_convergencia_log(historicos, f_otimo)
% Grafico proprio, em ESCALA LOGARITMICA. O plot_convergencia compartilhado
% usa eixo linear, inutil aqui: f cai de ~1e-2 para ~1e-12 e todas as curvas
% colapsariam sobre o zero.

fig = figure('Name', 'Convergencia — Trem de Engrenagens', 'Color', 'w');
hold on; grid on;
for r = 1:numel(historicos)
    h = historicos{r};
    h(~isfinite(h) | h <= 0) = NaN;    % Inf das particulas fora do dominio
    semilogy(h, 'Color', [0.6 0.6 0.85], 'LineWidth', 0.8);
end
yline_val = f_otimo;
plot([1 max(cellfun(@numel, historicos))], [yline_val yline_val], ...
     'r--', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
xlabel('Iteracao'); ylabel('f(x) = [1/6.931 - za zb/(zc zd)]^2');
title('PSO-RID — Trem de Engrenagens [DF2011] Sec. 5.1');
legend({'execucoes', 'otimo global'}, 'Location', 'northeast');
hold off;
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
