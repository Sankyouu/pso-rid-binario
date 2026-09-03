function resultados = run_todos_testes()
% RUN_TODOS_TESTES  Executa toda a suite de testes unitarios do projeto CPIO III.
%
% Framework: matlab.unittest (nativo do MATLAB).
%
% ATENCAO DE AMBIENTE: o GNU Octave NAO implementa matlab.unittest. Esta
% suite exige MATLAB. O Octave continua servindo para rodar os solvers
% (pso_rid, fem_*) diretamente.
%
% USO:
%   run_todos_testes                 % executa e imprime o resumo
%   r = run_todos_testes;            % devolve o array de TestResult
%
% Da linha de comando:
%   matlab -batch "cd('06_testes'); run_todos_testes"
%
% ORGANIZACAO DA SUITE (espelha os 3 blocos do projeto):
%   Bloco 1 - PSO-RID
%       TestDecodificadorPSO    Eq. (5), Eq. (6), Sec. 4.1  [D1][D2]
%       TestRegrasDeb           Deb (2000) pag. 316
%       TestVelocidadeBinaria   Eq. (7)-(8), Tabelas 1-2    [D3][D4]
%       TestMutacaoPolinomial   Eq. (9)                     [D5]
%   Bloco 2 - FEM
%       TestFemLinear           solucoes fechadas, simetria, equilibrio
%       TestFemNaoLinear        validacao analitica exata, K_T = K_E + K_G
%   Bloco 1+2+3
%       TestIntegracao          interfaces entre blocos, ciclo completo
%       TestOrquestradoresNovos reproducao dos valores publicados nos 4
%                               casos novos (20 e 51 barras, engrenagens, mola)
%
% NOTA: os casos de estudo nao sao arquivos proprios — vivem como funcao
% local no respectivo orquestrador e sao obtidos por acessor
% (main_hadi_nao_linear('caso'), main_awruch_discreto('caso')). Por isso
% 03_orquestrador precisa estar no path para os testes do Bloco 2.
% O oraculo analitico problema_trelica_rasa_2barras.m fica nesta pasta.

pasta_testes = fileparts(mfilename('fullpath'));
raiz         = fullfile(pasta_testes, '..');

addpath(genpath(fullfile(raiz, '01_pso_rid')));
addpath(genpath(fullfile(raiz, '02_fem_nao_linear')));
addpath(genpath(fullfile(raiz, '03_orquestrador')));
addpath(fullfile(raiz, '05_legado'));
addpath(pasta_testes);
addpath(raiz);

fprintf('==========================================================\n');
fprintf(' SUITE DE TESTES — CPIO III\n');
fprintf(' Pasta: %s\n', pasta_testes);
fprintf('==========================================================\n\n');

resultados = runtests(pasta_testes);

n_total  = numel(resultados);
n_ok     = sum([resultados.Passed]);
n_falhou = sum([resultados.Failed]);
n_incomp = sum([resultados.Incomplete]);
duracao  = sum([resultados.Duration]);

fprintf('\n==========================================================\n');
fprintf(' RESUMO\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Total      : %d\n', n_total);
fprintf(' Passaram   : %d\n', n_ok);
fprintf(' Falharam   : %d\n', n_falhou);
fprintf(' Incompletos: %d\n', n_incomp);
fprintf(' Duracao    : %.2f s\n', duracao);
fprintf('==========================================================\n');

if n_falhou > 0
    fprintf('\nTESTES QUE FALHARAM:\n');
    for k = find([resultados.Failed])
        fprintf('  - %s\n', resultados(k).Name);
    end
end

fprintf('\nFIM_TESTES\n');

end
