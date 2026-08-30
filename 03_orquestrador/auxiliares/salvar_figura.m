function salvar_figura(fig_handle, nome_arquivo)
% SALVAR_FIGURA Salva figuras nas pastas de resultados em PNG e FIG
%
% Exemplo:
%   salvar_figura(gcf, 'convergencia_hadi_100_particulas')

caminho_base = fullfile(fileparts(mfilename('fullpath')), '..', '..', '04_resultados', 'figuras');
if ~exist(caminho_base, 'dir')
    mkdir(caminho_base);
end

caminho_png = fullfile(caminho_base, [nome_arquivo, '.png']);
caminho_fig = fullfile(caminho_base, [nome_arquivo, '.fig']);

saveas(fig_handle, caminho_png);
saveas(fig_handle, caminho_fig);

fprintf('Figura salva com sucesso em:\n -> %s\n -> %s\n', caminho_png, caminho_fig);

end
