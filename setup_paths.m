function setup_paths(verbose)
%% SETUP_PATHS Adiciona as pastas do projeto CPIO III ao PATH do MATLAB/Octave
%
% ORGANIZACAO EM 3 BLOCOS
% (ver 00_docs/notas_e_relatorios/plano_final_reorganizacao_2026-08-28.md)
%
%   Bloco 1 — PSO-RID .......... 01_pso_rid/
%       pso_rid.m e AUTOCONTIDO: o loop e todas as pecas do algoritmo (uma
%       funcao local por formula do artigo) estao no mesmo arquivo.
%       Handles das funcoes locais: aux = pso_rid('auxiliares').
%
%   Bloco 2 — FEM Nao Linear ... 02_fem_nao_linear/
%       SO os solvers genericos, sem dados de nenhum problema:
%       fem_linear_solver.m e fem_nao_linear_solver.m.
%
%   Bloco 3 — Orquestrador ..... 03_orquestrador/
%       main_*.m       : um experimento por arquivo, autocontido, com a
%                        DEFINICAO DO PROBLEMA embutida como funcao local
%       auxiliares/    : apenas blocos grandes de impressao/grafico
%
%       Casos de estudo (acessores, sem rodar otimizacao):
%           caso = main_hadi_nao_linear('caso')   benchmark [HA2003]
%           caso = main_awruch_discreto('caso')   catalogos por barra
%
%   Testes ..................... 06_testes/  (matlab.unittest — exige MATLAB)
%       Inclui problema_trelica_rasa_2barras.m, o oraculo analitico usado
%       apenas pelos testes do solver (nao e um caso de estudo).
%   Legado ..................... 05_legado/  (nao entra no path)
%
% Uso:
%   setup_paths          % adiciona os caminhos e exibe as pastas
%   setup_paths(false)   % adiciona silenciosamente

if nargin < 1
    verbose = true;
end

root_dir = fileparts(mfilename('fullpath'));

dirs = { ...
    fullfile(root_dir, '01_pso_rid')        % Bloco 1
    fullfile(root_dir, '02_fem_nao_linear') % Bloco 2 (solver + problemas)
    fullfile(root_dir, '03_orquestrador')   % Bloco 3 (+ auxiliares)
    };

for k = 1:numel(dirs)
    if exist(dirs{k}, 'dir')
        addpath(genpath(dirs{k}));
    end
end

if verbose
    fprintf('[OK] Caminhos do projeto CPIO III adicionados ao path:\n');
    for k = 1:numel(dirs)
        if exist(dirs{k}, 'dir')
            fprintf('  -> %s\n', dirs{k});
        end
    end
    fprintf('\nExperimentos disponiveis em 03_orquestrador/:\n');
    fprintf('  main_hadi_nao_linear      Benchmark Hadi 2003, analise nao linear\n');
    fprintf('  main_hadi_linear          Benchmark Hadi 2003, analise linear\n');
    fprintf('  main_awruch_discreto      Catalogos independentes por barra\n');
    fprintf('  main_estudo_estatistico   Comparacao multi-semente de configuracoes\n');
    fprintf('\nTestes: matlab -batch "cd(''06_testes''); run_todos_testes"\n');
end
end
