function garantir_caminhos(nomes_requeridos)
% GARANTIR_CAMINHOS  Adiciona as pastas do projeto ao path, se necessario.
%
% Usado no topo de cada orquestrador (03_orquestrador/main_*.m) para que
% funcione mesmo se chamado diretamente, sem setup_paths ter rodado antes.
% Centraliza o que antes era uma funcao local quase identica, copiada em
% cada um dos 7 orquestradores.
%
% ENTRADA
%   nomes_requeridos : cell array de nomes de função a checar via
%                      exist(nome,'file')==2 (ex.: {'pso_rid','fem_nao_linear_solver'})
%
% See also setup_paths

if all(cellfun(@(n) exist(n, 'file') == 2, nomes_requeridos))
    return;
end

% Este arquivo vive em 03_orquestrador/auxiliares/, um nivel ABAIXO dos
% orquestradores — por isso '..','..' ate a raiz do projeto. O padrao
% antigo (copiado dentro de 03_orquestrador/) usava so '..', porque cada
% copia vivia um nivel acima daqui.
raiz = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(raiz);
setup_paths(false);

end
