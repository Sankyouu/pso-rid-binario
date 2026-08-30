function [mapa, total_dim] = rid_mapear_dimensoes(config_vars)
% RID_MAPEAR_DIMENSOES  Define quantas dimensoes cada variavel ocupa na particula.
%
% BLOCO 1 (PSO-RID) — auxiliar isolado e testavel.
%
% REFERENCIA:
%   [DF2011] Datta, D.; Figueira, J.R. "A real-integer-discrete-coded particle
%            swarm optimization for design problems".
%            Applied Soft Computing 11 (2011) 3625-3633.
%
% Eq. (5) [DF2011]:   L = R + sum_j B_j^(I) + sum_k B_k^(D)
%   L       = numero total de dimensoes posicionais da particula
%   R       = numero de variaveis reais (1 dimensao cada)
%   B_j^(I) = bits da j-esima variavel inteira
%   B_k^(D) = bits da k-esima variavel discreta
%
% [DF2011] Sec. 4.1:
%   "The number of dimensions required for an integer variable is determined
%    from the maximum number of binary bits required to represent its upper
%    limit. It is to be mentioned that a discrete variable is also dealt with
%    as an integer variable, an integer value of which represents the index of
%    its actual discrete value, so that the lower limit of such a variable is 1
%    and the upper limit is the number of its allowable discrete values."
%
% DIMENSIONAMENTO DE BITS (deduzido da Eq. 6 e conferido com os exemplos
% numericos do proprio artigo):
%
%   INTEIRO  : a Eq. (6) precisa gerar o proprio valor var.max, logo
%              2^B - 1 >= max  =>  B = ceil(log2(max + 1)).
%              Confere com Sec. 5.1 (trem de engrenagens): 12 <= z <= 60
%              com "6 binary bits" -> ceil(log2(61)) = 6.
%              Confere com Sec. 5.2 (mola): N inteiro em [1,70] com
%              "7 binary bits" -> ceil(log2(71)) = 7.
%
%   DISCRETO : o indice vai de 1 a N, ou seja x_int (Eq. 6) vai de 0 a N-1,
%              logo 2^B >= N  =>  B = ceil(log2(N)).
%              Confere com Sec. 5.2 (mola): 42 diametros discretos com
%              "6 binary bits" -> ceil(log2(42)) = 6.
%
% ENTRADA
%   config_vars : struct array com .tipo ('R','I','D') e os campos exigidos
%                 por cada tipo (.min/.max para R e I; .opcoes para D)
%
% SAIDAS
%   mapa      : struct array com, para cada variavel,
%                 .idx_original, .tipo, .n_bits, .inicio, .fim
%   total_dim : numero total de dimensoes da particula (L da Eq. 5)

template.idx_original = 0;
template.tipo         = '';
template.n_bits       = 0;
template.inicio       = 0;
template.fim          = 0;

mapa      = repmat(template, length(config_vars), 1);
total_dim = 0;

for i = 1:length(config_vars)
    var = config_vars(i);

    info.idx_original = i;
    info.tipo         = var.tipo;

    if strcmp(var.tipo, 'R')
        % Variavel real: 1 dimensao, sem bits. Eq. (5) [DF2011], termo R.
        info.n_bits = 0;
        info.inicio = total_dim + 1;
        total_dim   = total_dim + 1;
        info.fim    = total_dim;
    else
        if strcmp(var.tipo, 'I')
            info.n_bits = max(ceil(log2(var.max + 1)), 1);
        else  % 'D'
            info.n_bits = max(ceil(log2(length(var.opcoes))), 1);
        end
        info.inicio = total_dim + 1;
        total_dim   = total_dim + info.n_bits;
        info.fim    = total_dim;
    end

    mapa(i) = info;
end

end
