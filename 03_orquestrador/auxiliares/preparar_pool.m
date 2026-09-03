function n = preparar_pool(n_pedido, n_tarefas)
% PREPARAR_POOL  Garante um pool paralelo com no maximo n_pedido workers,
% respeitando um teto por memoria disponivel. Devolve quantos workers usar
% de fato no parfor (0 = execucao serial, sem pool nenhum).
%
% Padrao extraido de 03_orquestrador/main_estudo_estatistico.m (funcoes
% locais preparar_pool/resolver_workers), reaproveitado aqui como helper
% compartilhado por 01_pso_rid/pso_rid.m (avaliacao de particulas em
% paralelo) e pelos orquestradores (laco multi-start em paralelo).
% main_estudo_estatistico.m NAO foi editado — so serviu de referencia.
%
% -------------------------------------------------------------------------
% POR QUE O TETO E POR MEMORIA, E NAO POR feature('numcores')
% -------------------------------------------------------------------------
% Um worker do Parallel Computing Toolbox nao e uma thread: e um processo
% MATLAB inteiro, com seu proprio interpretador e copia dos dados fatiados.
% Em regime, cada um ocupa da ordem de 2 GB residentes. Numa maquina de
% muitos nucleos mas pouca RAM livre, um pool "um worker por nucleo" pede
% memoria demais e trava o sistema, sobretudo sem swap. O limite certo e
% quantos processos cabem na memoria que sobra, nao o numero de nucleos.
%
% -------------------------------------------------------------------------
% ENTRADAS
% -------------------------------------------------------------------------
%   n_pedido  : numero de workers pedido pelo chamador. 0 (ou negativo)
%               forca execucao serial, sem abrir pool nenhum.
%   n_tarefas : (opcional) numero de tarefas do laco que vai usar o pool —
%               pedir mais workers do que tarefas so desperdica memoria.
%               Default: Inf (sem limite adicional).
%
% SAIDA
%   n : quantos workers usar de fato em parfor(i=1:N, n). 0 = serial.

if nargin < 2 || isempty(n_tarefas), n_tarefas = inf; end

if n_pedido <= 0
    n = 0;
    return;
end

if isempty(ver('parallel')) || ~license('test', 'Distrib_Computing_Toolbox')
    n = 0;
    return;
end

% 2 GB/worker de proposito, com folga: errar para cima aqui custa alguns
% workers a menos; errar para baixo trava a maquina. RESERVA_GB cobre o SO e
% a sessao MATLAB cliente. TETO_WORKERS e um teto rigido — acima disso o
% ganho de paralelismo nao paga o risco de memoria.
RAM_POR_WORKER_GB = 2;
RESERVA_GB        = 3;
TETO_WORKERS      = 6;

livre_gb = memoria_disponivel_gb();
if isnan(livre_gb)
    % Nao deu para medir a memoria: escolhe o minimo util, nao o maximo
    % possivel — 2 workers cabem em praticamente qualquer maquina que rode
    % MATLAB, e o chamador continua CORRETO (so mais lento) se sobrar pouco.
    n_por_memoria = 2;
else
    n_por_memoria = floor(max(livre_gb - RESERVA_GB, 0) / RAM_POR_WORKER_GB);
end

n = min([n_pedido, n_por_memoria, max(feature('numcores') - 1, 1), ...
         TETO_WORKERS, n_tarefas]);

% Abaixo de 2 workers o paralelismo so custa: pagar-se-ia a abertura do pool
% e a serializacao das tarefas para rodar praticamente em serie.
if n < 2
    n = 0;
    return;
end

pool = gcp('nocreate');

if isempty(pool)
    fprintf('    Abrindo pool com %d workers (~%d GB estimados)...\n', n, 2*n);
    parpool('Processes', n);
    return;
end

% Ja havia um pool: reaproveita, sem derrubar nem redimensionar — isso
% mataria um pool que o usuario possa estar usando para outra coisa.
if pool.NumWorkers > n
    fprintf(['    Pool existente tem %d workers (pedido: %d). ' ...
             'Reaproveitando; o parfor usa no maximo %d.\n'], ...
            pool.NumWorkers, n, n);
else
    n = pool.NumWorkers;
    fprintf('    Reaproveitando pool existente com %d workers.\n', n);
end

end


function gb = memoria_disponivel_gb()
% MEMORIA_DISPONIVEL_GB  Memoria utilizavel agora, em GB. NaN se nao der
% para medir.
%
% MemAvailable do /proc/meminfo e a estimativa do proprio kernel de quanto
% pode ser alocado sem entrar em swap — ja inclui o cache que ele devolveria
% sob pressao. E a medida certa aqui; MemFree subestimaria muito.

gb = NaN;

try
    if isunix && ~ismac
        txt = fileread('/proc/meminfo');
        tok = regexp(txt, 'MemAvailable:\s+(\d+)\s+kB', 'tokens', 'once');
        if ~isempty(tok)
            gb = str2double(tok{1}) / 1024 / 1024;
        end
    elseif ispc
        [~, sys] = memory();
        gb = sys.PhysicalMemory.Available / 1024^3;
    end
catch
    gb = NaN;   % qualquer falha cai no caminho conservador de quem chamou
end
end
