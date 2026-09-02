function resultado = main_estudo_estatistico(n_sementes, configuracoes, tipo_analise)
% MAIN_ESTUDO_ESTATISTICO  Comparacao multi-semente de configuracoes do PSO-RID.
%
% -------------------------------------------------------------------------
% MOTIVACAO
% -------------------------------------------------------------------------

% Este orquestrador compara configuracoes do PSO-RID sobre o MESMO conjunto
% de sementes (amostras pareadas), reportando media, desvio padrao,
% melhor/pior e taxa de viabilidade. Foi criado para responder duas questoes
% que execucoes unicas nao resolviam, dada a variancia do algoritmo.
%
% -------------------------------------------------------------------------
% QUESTOES JA RESPONDIDAS POR ESTE ESTUDO
% -------------------------------------------------------------------------

%   Q1. A decodificacao discreta literal de [DF2011] ('datta') e melhor ou
%       pior que a heuristica 'proporcional'?
%       RESPONDIDA (30 sementes pareadas, FEM linear, 2026-08-31):
%         'datta'        melhor 2357.97 | media 2572.64 | desvio 173.05
%         'proporcional' melhor 2342.11 | media 2504.72 | desvio  97.40
%       Variancia: F = 3.16 (critico 2.10) -> 'proporcional'
%       SIGNIFICATIVAMENTE mais consistente. Media: t = -1.836 (sem
%       evidencia), mas teste do sinal 21/30 (p = 0.043, significativo).
%       DECISAO: 'proporcional' virou o padrao — ver [D1] em pso_rid.m.
%
%   Q2. A auto-adaptacao instantanea dos hiperparametros ([DF2011] Sec. 5)
%       prejudica a convergencia frente a valores fixos?
%       RESPONDIDA (10 sementes pareadas): +31.30 kg, t = 0.631, 4/10
%       vitorias. SEM evidencia. A suspeita original vinha de uma unica
%       execucao (legado 2356.80 vs nova 2536.56) e era artefato de
%       amostra de tamanho 1. HIPOTESE REFUTADA.
%
% As configuracoes padrao abaixo foram reorientadas para VIGIAR essas
% decisoes, nao mais para descobri-las.
%
% -------------------------------------------------------------------------
% METODOLOGIA
% -------------------------------------------------------------------------
%
% - Cada configuracao e executada uma vez por semente.
% - Todas as configuracoes veem EXATAMENTE as mesmas sementes (pareamento),
%   o que reduz a variancia da comparacao.
% - O orcamento e controlado por numero de AVALIACOES FEM, nao por
%   iteracoes: configuracoes com muito descarte estrutural (modo 'datta')
%   fariam mais iteracoes no mesmo tempo, o que enviesaria a comparacao por
%   iteracao. [DF2011] Sec. 5 tambem mede custo em numero de avaliacoes da
%   funcao objetivo.
% - Reporta-se o melhor peso VIAVEL de cada execucao (violacao <= 0).
%
% ATENCAO AO CUSTO COMPUTACIONAL: cada avaliacao chama o solver nao linear.
% Com os padroes abaixo (4 configuracoes x 10 sementes) o estudo leva horas.
% Use n_sementes menor para um ensaio rapido.
%
% -------------------------------------------------------------------------
% USO
% -------------------------------------------------------------------------
%
%   main_estudo_estatistico                  % 10 sementes, 4 configuracoes
%   main_estudo_estatistico(3)               % ensaio rapido
%   main_estudo_estatistico(10, [], 'linear')% FEM linear (muito mais rapido)
%   r = main_estudo_estatistico(5);
%
% ENTRADAS
%
%   n_sementes    : quantas execucoes por configuracao      (padrao 10)
%   configuracoes : struct array com .nome e .params        (padrao: as 4 abaixo)
%   tipo_analise  : 'nao_linear' (padrao) ou 'linear'
%
% SAIDA
%
%   resultado : struct array com as estatisticas de cada configuracao

if nargin < 1 || isempty(n_sementes),   n_sementes   = 10;           end
if nargin < 3 || isempty(tipo_analise), tipo_analise = 'nao_linear'; end

garantir_caminhos();

% Caso Hadi — fonte unica: funcao local de main_hadi_nao_linear.m,
% exposta pelo acessor (ver cabecalho daquele arquivo).

caso = main_hadi_nao_linear('caso');

% -------------------------------------------------------------------------
% CONFIGURACOES COMPARADAS
% -------------------------------------------------------------------------
% Desenhadas para isolar UMA variavel por vez em relacao a configuracao
% base (fiel a [DF2011]):
%   A: PADRAO ATUAL do solver ('proporcional' + auto-adaptativo)
%   B: A, trocando so a decodificacao para 'datta'  -> reexamina [D1]
%   C: A, trocando so a auto-adaptacao por fixos    -> reexamina Q2
%
% NOTA DE HISTORICO (2026-08-31): ate esta data a base era 'datta', e o
% estudo existia para responder Q1 (decodificacao) e Q2 (auto-adaptacao).
% Ambas foram respondidas — ver o relatorio em 00_docs/notas_e_relatorios/ e
% o bloco [D1] em pso_rid.m. A base passou a ser 'proporcional', e as
% configuracoes abaixo foram reorientadas para VIGIAR essas decisoes: se
% alguma delas deixar de se sustentar num novo desenho experimental, aparece
% aqui.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(configuracoes)
    base                           = struct();
    base.n_particulas              = 100;
    base.w                         = 1.0;
    base.c1                        = 1.0;
    base.c2                        = 2.0;
    base.pm                        = 0.15;
    base.auto_adaptativo           = true;
    base.decodificacao_discreta    = 'proporcional';   % [D1]
    base.tol_estagnacao            = 200;
    base.verbose                   = false;

    cfgA = base;

    cfgB = base;
    cfgB.decodificacao_discreta    = 'datta';

    cfgC = base;
    cfgC.auto_adaptativo           = false;
    cfgC.w                         = 0.9;
    cfgC.c1                        = 2.0;
    cfgC.c2                        = 2.0;

    configuracoes = struct( ...
        'nome',   {'A: padrao atual (proporcional)', ...
                   'B: A + decod. datta', ...
                   'C: A + hiperparam. fixos'}, ...
        'params', {cfgA, cfgB, cfgC});
end

% Orcamento em avaliacoes FEM (igual para todas as configuracoes)
orcamento_avaliacoes = 15000;

sementes = 1:n_sementes;

imprimir_cabecalho(caso, configuracoes, sementes, orcamento_avaliacoes, tipo_analise);

% -------------------------------------------------------------------------
% EXECUCAO
% -------------------------------------------------------------------------

n_cfg     = numel(configuracoes);
pesos     = nan(n_cfg, n_sementes);
viaveis   = false(n_cfg, n_sementes);
avaliacoes= nan(n_cfg, n_sementes);
tempos    = nan(n_cfg, n_sementes);

for c = 1:n_cfg
    fprintf('\n>>> Configuracao %d/%d — %s\n', c, n_cfg, configuracoes(c).nome);

    for s = 1:n_sementes
        params = configuracoes(c).params;
        % max_iter alto: a parada efetiva e o orcamento de avaliacoes,
        % imposto pelo envoltorio contador abaixo.
        params.max_iter = 100000;

        rng(sementes(s));

        % Funcao objetivo com orcamento de avaliacoes do FEM. O contador e o
        % melhor-ate-agora vivem numa closure (ver criar_funcao_com_orcamento).
        [fobj, ler_melhor] = criar_funcao_com_orcamento( ...
                                 caso, tipo_analise, orcamento_avaliacoes);

        t0 = tic;
        try
            % Caminho normal: o PSO termina por estagnacao ou max_iter antes
            % de esgotar o orcamento.
            [~, peso_s, ~, det_s] = pso_rid(fobj, caso.config_vars, params);
            viol_s = det_s.gbest_viol;
            [~, ~, ~, n_eval] = ler_melhor();
        catch err
            % Caminho esperado quando o orcamento acaba primeiro: recupera o
            % melhor resultado registrado pela closure antes da interrupcao.
            if ~strcmp(err.identifier, 'main_estudo_estatistico:orcamentoEsgotado')
                rethrow(err);
            end
            [~, peso_s, viol_s, n_eval] = ler_melhor();
        end
        tempos(c,s) = toc(t0);

        pesos(c,s)      = peso_s;
        viaveis(c,s)    = (viol_s <= 0);
        avaliacoes(c,s) = n_eval;

        fprintf('    semente %2d: peso = %8.2f kg | viavel = %d | avaliacoes = %5d | %.0fs\n', ...
                sementes(s), peso_s, viaveis(c,s), avaliacoes(c,s), tempos(c,s));
    end
end

% -------------------------------------------------------------------------
% ESTATISTICAS E RELATORIO
% -------------------------------------------------------------------------

resultado = struct([]);
for c = 1:n_cfg
    p_viaveis = pesos(c, viaveis(c,:));

    resultado(c).nome            = configuracoes(c).nome;
    resultado(c).params          = configuracoes(c).params;
    resultado(c).pesos           = pesos(c,:);
    resultado(c).viaveis         = viaveis(c,:);
    resultado(c).taxa_viabilidade= mean(viaveis(c,:));
    resultado(c).avaliacoes      = avaliacoes(c,:);
    resultado(c).tempo_medio     = mean(tempos(c,:));

    if isempty(p_viaveis)
        [resultado(c).melhor, resultado(c).media, ...
         resultado(c).desvio, resultado(c).pior] = deal(NaN);
    else
        resultado(c).melhor = min(p_viaveis);
        resultado(c).media  = mean(p_viaveis);
        resultado(c).desvio = std(p_viaveis);
        resultado(c).pior   = max(p_viaveis);
    end
end

imprimir_tabela_comparativa(resultado, caso.ref_peso, n_sementes);
imprimir_conclusoes(resultado);
salvar_log_estudo(resultado, sementes, orcamento_avaliacoes, tipo_analise);

end


% #########################################################################
% FUNCOES LOCAIS
% #########################################################################

function [fobj, ler_melhor] = criar_funcao_com_orcamento(caso, tipo_analise, limite)

% Cria uma funcao objetivo com ORCAMENTO DE AVALIACOES do FEM, para que
% configuracoes com taxas de descarte estrutural diferentes sejam comparadas
% de forma justa (ver METODOLOGIA no topo do arquivo).
%
% MECANISMO DE PARADA
% -------------------
% Ao esgotar o orcamento, a funcao LANCA um erro identificado
% ('main_estudo_estatistico:orcamentoEsgotado'), que o chamador captura para
% encerrar a execucao imediatamente.
%
% Por que lancar erro em vez de devolver custo infinito: devolvendo Inf, o
% PSO continuaria iterando ate max_iter (ou ate a estagnacao), gastando
% tempo sem realizar nenhuma avaliacao util. Depender da estagnacao para
% terminar seria implicito e fragil — bastaria alguem afrouxar
% tol_estagnacao para o laco girar dezenas de milhares de iteracoes a toa.
%
% Como o erro interrompe o pso_rid antes do seu retorno normal, esta closure
% mantem por conta propria o MELHOR RESULTADO ja visto, aplicando a mesma
% regra de dominancia de Deb usada internamente pelo solver (domina_deb).
% Assim nada se perde ao abortar.

% domina_deb e funcao LOCAL de pso_rid.m; o handle vem pelo acessor
% documentado no cabecalho daquele arquivo. Usar exatamente a mesma
% implementacao (e nao uma copia) e o que garante que o "melhor" registrado
% aqui coincida com o criterio interno do otimizador.

aux_pso     = pso_rid('auxiliares');
domina_deb  = aux_pso.domina_deb;

n_chamadas    = 0;
melhor_areas  = [];
melhor_custo  = inf;
melhor_viol   = inf;

    function [custo, violacao] = avaliar(areas)
        if n_chamadas >= limite
            error('main_estudo_estatistico:orcamentoEsgotado', ...
                  'Orcamento de %d avaliacoes do FEM esgotado.', limite);
        end
        n_chamadas = n_chamadas + 1;

        [custo, violacao] = avaliar_projeto(areas, caso, tipo_analise);

        % Registra o melhor pela regra de [DEB2000] (mesma do solver)
        if domina_deb(custo, violacao, melhor_custo, melhor_viol)
            melhor_custo = custo;
            melhor_viol  = violacao;
            melhor_areas = areas;
        end
    end

    function [areas, custo, violacao, n] = ler()
        areas    = melhor_areas;
        custo    = melhor_custo;
        violacao = melhor_viol;
        n        = n_chamadas;
    end

fobj       = @avaliar;
ler_melhor = @ler;
end


function [custo, violacao] = avaliar_projeto(areas, caso, tipo_analise)
% Funcao objetivo — mesma formulacao dos demais orquestradores.
% [HA2003] Eq. (1) objetivo, Eq. (2) restricoes. Violacao devolvida
% separadamente para uso com a regra de [DEB2000].

if strcmp(tipo_analise, 'linear')
    [peso, Sigma, u_livre] = fem_linear_solver(caso, areas);
else
    [peso, Sigma, u_livre] = fem_nao_linear_solver(caso, areas);
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


function imprimir_cabecalho(caso, configuracoes, sementes, orcamento, tipo)
fprintf('\n================================================================\n');
fprintf(' ESTUDO ESTATISTICO MULTI-SEMENTE — PSO-RID\n');
fprintf(' Problema: %s\n', caso.nome);
fprintf('----------------------------------------------------------------\n');
fprintf(' Analise estrutural       : %s\n', tipo);
fprintf(' Configuracoes comparadas : %d\n', numel(configuracoes));
fprintf(' Sementes por configuracao: %d  (%s)\n', numel(sementes), mat2str(sementes));
fprintf(' Orcamento por execucao   : %d avaliacoes FEM\n', orcamento);
fprintf(' Referencia [HA2003]      : %.2f kg\n', caso.ref_peso);
fprintf('----------------------------------------------------------------\n');
for c = 1:numel(configuracoes)
    p = configuracoes(c).params;
    fprintf('  %s\n', configuracoes(c).nome);
    fprintf('     decodificacao=%s | auto_adapt=%d | w=%.2f c1=%.2f c2=%.2f\n', ...
            p.decodificacao_discreta, p.auto_adaptativo, p.w, p.c1, p.c2);
end
fprintf('================================================================\n');
end


function imprimir_tabela_comparativa(res, ref_peso, n_sementes)
fprintf('\n================================================================\n');
fprintf(' RESULTADOS (peso em kg, apenas execucoes VIAVEIS)\n');
fprintf('----------------------------------------------------------------\n');
fprintf(' %-34s %8s %8s %8s %6s\n', 'Configuracao', 'Melhor', 'Media', 'Desvio', 'Viav.');
fprintf(' %s\n', repmat('-', 1, 68));
for c = 1:numel(res)
    fprintf(' %-34s %8.2f %8.2f %8.2f %5.0f%%\n', ...
            res(c).nome, res(c).melhor, res(c).media, res(c).desvio, ...
            100*res(c).taxa_viabilidade);
end
fprintf(' %s\n', repmat('-', 1, 68));
if ~isnan(ref_peso)
    fprintf(' Referencia [HA2003] Tab. 1 (Case 2 discreto): %.2f kg\n', ref_peso);
end
fprintf(' Amostras por configuracao: %d sementes pareadas\n', n_sementes);
fprintf('================================================================\n');
end


function imprimir_conclusoes(res)

% Comparacoes PAREADAS que respondem Q1 e Q2 da documentacao no topo.
%
% METODOLOGIA — por que analisar as DIFERENCAS, e nao as medias:
% as configuracoes sao executadas sobre as MESMAS sementes, formando
% amostras pareadas. Nesse desenho, comparar a diferenca das medias contra
% os desvios MARGINAIS de cada configuracao e o teste errado: os desvios
% marginais carregam a variabilidade entre sementes, que o pareamento
% justamente elimina.
%
% O correto e trabalhar com o vetor de diferencas por semente
%     d_i = peso_config2(i) - peso_config1(i)
% e comparar media(d) com o ERRO PADRAO da media, desvio(d)/sqrt(n).
% Isso e a estatistica t pareada. Reporta-se tambem a contagem de vitorias,
% que e robusta a valores extremos.

fprintf('\n LEITURA DAS COMPARACOES PAREADAS\n');
fprintf('----------------------------------------------------------------\n');

% Os rotulos supoem as configuracoes PADRAO. Com configuracoes proprias,
% cada bloco compara res(k) contra res(1) — leia pelos nomes impressos.
rotulos = {'[D1] decodificacao: datta vs padrao proporcional', ...
           'Q2 — hiperparametros fixos vs auto-adaptativos'};
for k = 2:min(numel(res), 3)
    comparar_pareado(res(k), res(1), rotulos{k-1});
end

fprintf('----------------------------------------------------------------\n');
fprintf(' NOTA: a estatistica t pareada acima pressupoe diferencas\n');
fprintf(' aproximadamente normais. Para um teste livre de distribuicao,\n');
fprintf(' aplique Wilcoxon pareado (signrank) sobre resultado(c).pesos.\n');
fprintf('================================================================\n\n');
end


function comparar_pareado(res_novo, res_ref, titulo)
% Compara duas configuracoes usando as diferencas semente a semente.

d = res_novo.pesos(:) - res_ref.pesos(:);
d = d(isfinite(d));
n = numel(d);

fprintf('\n %s\n', titulo);
fprintf('   (%s)  -  (%s)\n', res_novo.nome, res_ref.nome);

if n < 2
    fprintf('   Amostras insuficientes para comparacao pareada.\n');
    return;
end

media_d  = mean(d);
desvio_d = std(d);
erro_pad = desvio_d / sqrt(n);

vitorias = sum(d < 0);   % configuracao nova produziu peso MENOR

fprintf('   media das diferencas : %+8.2f kg\n', media_d);
fprintf('   desvio das difs      : %8.2f kg\n', desvio_d);
fprintf('   erro padrao da media : %8.2f kg\n', erro_pad);
fprintf('   vitorias             : %d de %d sementes\n', vitorias, n);

if erro_pad > 0
    t = media_d / erro_pad;
    % Valor critico bilateral aproximado de t para 5%, gl = n-1
    t_crit = valor_critico_t(n - 1);
    fprintf('   t pareado            : %+8.3f  (critico ~%.2f a 5%%)\n', t, t_crit);

    if abs(t) < t_crit
        fprintf('   >> SEM evidencia estatistica de diferenca (|t| < critico).\n');
        fprintf('      A diferenca observada e compativel com ruido do algoritmo.\n');
    elseif media_d < 0
        fprintf('   >> Diferenca SIGNIFICATIVA: a configuracao nova produz\n');
        fprintf('      pesos menores.\n');
    else
        fprintf('   >> Diferenca SIGNIFICATIVA: a configuracao de referencia\n');
        fprintf('      produz pesos menores.\n');
    end
end

% A dispersao tambem importa na pratica: uma configuracao com media
% equivalente mas menor variabilidade e preferivel por ser mais confiavel.
if res_novo.desvio < res_ref.desvio * 0.7
    fprintf('   >> Observacao: a configuracao nova e MAIS CONSISTENTE\n');
    fprintf('      (desvio %.2f vs %.2f).\n', res_novo.desvio, res_ref.desvio);
elseif res_ref.desvio < res_novo.desvio * 0.7
    fprintf('   >> Observacao: a configuracao de referencia e MAIS CONSISTENTE\n');
    fprintf('      (desvio %.2f vs %.2f).\n', res_ref.desvio, res_novo.desvio);
end
end


function t_crit = valor_critico_t(gl)
% Valor critico bilateral de t a 5%, tabelado para os graus de liberdade
% mais usuais. Evita depender do Statistics Toolbox (tinv).
tabela_gl = [1 2 3 4 5 6 7 8 9 10 12 15 20 25 30 40 60 120];
tabela_t  = [12.71 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 ...
             2.228 2.179 2.131 2.086 2.060 2.042 2.021 2.000 1.980];

if gl <= 0
    t_crit = Inf;
elseif gl >= tabela_gl(end)
    t_crit = 1.96;                       % limite normal
else
    t_crit = interp1(tabela_gl, tabela_t, gl, 'linear');
end
end


function salvar_log_estudo(res, sementes, orcamento, tipo)
% Registra o estudo em 04_resultados/logs/ para consulta posterior.
pasta = fullfile(fileparts(mfilename('fullpath')), '..', '04_resultados', 'logs');
if ~exist(pasta, 'dir'); mkdir(pasta); end

nome = fullfile(pasta, sprintf('estudo_estatistico_%s_%s.txt', ...
                               tipo, datestr(now, 'yyyymmdd_HHMMSS')));

fid = fopen(nome, 'w');
if fid < 0
    warning('main_estudo_estatistico:logNaoSalvo', ...
        'Nao foi possivel gravar o log em %s.', nome);
    return;
end

fprintf(fid, 'ESTUDO ESTATISTICO PSO-RID\n');
fprintf(fid, 'Data: %s\n', datestr(now));
fprintf(fid, 'Analise: %s | Sementes: %s | Orcamento: %d avaliacoes\n\n', ...
        tipo, mat2str(sementes), orcamento);

for c = 1:numel(res)
    fprintf(fid, '%s\n', res(c).nome);
    fprintf(fid, '  melhor=%.2f media=%.2f desvio=%.2f viabilidade=%.0f%%\n', ...
            res(c).melhor, res(c).media, res(c).desvio, 100*res(c).taxa_viabilidade);
    fprintf(fid, '  pesos por semente: %s\n\n', mat2str(res(c).pesos, 6));
end
fclose(fid);

fprintf(' Log do estudo salvo em:\n   %s\n\n', nome);
end


function garantir_caminhos()
if exist('pso_rid', 'file') == 2 && exist('fem_nao_linear_solver', 'file') == 2
    return;
end
raiz = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(raiz);
setup_paths(false);
end
