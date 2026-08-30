function caso = problema_trelica_rasa_2barras(a, h, E, A, dens)
% PROBLEMA_TRELICA_RASA_2BARRAS  Trelica rasa simetrica com SOLUCAO ANALITICA EXATA.
%
% BLOCO 2b do projeto CPIO III — ARQUIVO DE PARAMETROS + oraculo analitico.
%
% Este caso existe para VALIDAR o solver nao linear contra uma solucao
% fechada, calculavel a mao, independente de qualquer outro artigo ou
% software. E o teste mais forte disponivel para o Newton-Raphson, porque
% nao depende de reproduzir o resultado de terceiros.
%
% =========================================================================
% GEOMETRIA (trelica de von Mises / trelica rasa de duas barras)
% =========================================================================
%
%                        no 3 (apice)
%                          o  <-- carga P, vertical para baixo
%                        /   \
%                      /   |   \            altura inicial = h
%                    /     |     \
%          no 1    o-------+-------o    no 2
%                  ^<--a-->|<--a-->^
%                 apoio           apoio
%
%   Ambos os apoios sao fixos. Por simetria o apice so se desloca na
%   vertical, o que reduz o problema a UMA incognita escalar (w) e permite
%   solucao fechada.
%
% =========================================================================
% SOLUCAO ANALITICA EXATA  (derivada para DEFORMACAO DE ENGENHARIA)
% =========================================================================
% Seja w o deslocamento vertical do apice, positivo para BAIXO.
%
%   comprimento inicial :  L0   = sqrt(a^2 + h^2)
%   comprimento atual   :  L(w) = sqrt(a^2 + (h - w)^2)
%   deformacao          :  eps  = (L - L0) / L0            (engenharia)
%   forca axial         :  N    = E * A * eps              (tracao > 0)
%
% Equilibrio VERTICAL do no do apice na configuracao DEFORMADA:
%   cada barra puxa/empurra o apice com componente vertical N*(h-w)/L,
%   e ha duas barras simetricas, logo
%
%          P(w) = -2 * N(w) * (h - w) / L(w)                        (*)
%
%   com P positivo para baixo.
%
% VERIFICACOES DE SANIDADE DA FORMULA (*):
%   - Para 0 < w < h temos L < L0, logo N < 0 (compressao) e P > 0:
%     carga para baixo produz deslocamento para baixo, como esperado.
%   - Em w = h o apice atinge a linha dos apoios, (h-w) = 0 e P = 0:
%     este e o ponto de snap-through classico da trelica de von Mises.
%
% CONSEQUENCIA PARA O USO: a estrutura tem PONTO LIMITE. Um solver de
% controle de CARGA (como o do Bloco 2a) so consegue seguir o caminho de
% equilibrio ABAIXO da carga critica. Os testes devem usar cargas com folga
% em relacao ao maximo de P(w) — o campo caso.P_critica abaixo fornece esse
% valor para dimensionar os ensaios.
%
% NOTA SOBRE A LITERATURA: formulas classicas de trelica rasa (por exemplo
% em Crisfield Vol. 1, e reproduzidas por Yaw, "3D Co-rotational Truss
% Formulation", Eq. 10.1) sao derivadas com DEFORMACAO DE GREEN e com
% hipotese de abatimento (shallowness). Elas NAO coincidem com (*) fora do
% regime de deformacoes infinitesimais. Como o solver do Bloco 2a usa
% deformacao de ENGENHARIA, a formula (*) — derivada aqui para essa mesma
% medida — e o oraculo correto. Usar a formula de Green produziria um
% desacordo que seria erroneamente interpretado como bug no solver.
%
% =========================================================================
% ENTRADAS (todas opcionais; os padroes formam um caso raso bem-condicionado)
% =========================================================================
%   a    : meio-vao horizontal            (padrao 1000  [mm])
%   h    : altura inicial do apice        (padrao  100  [mm])
%   E    : modulo de elasticidade         (padrao 70000 [MPa])
%   A    : area da secao de cada barra    (padrao  100  [mm^2])
%   dens : densidade                      (padrao 2.7e-6 [kg/mm^3])
%
% SAIDA
%   caso : struct com os campos exigidos pelos solvers do Bloco 2a
%          (.nodes0 .elements .apoios .F_total .E .dens) mais:
%          .P_de_w      handle da solucao analitica (*)
%          .w_de_P      handle que inverte (*) numericamente
%          .P_critica   carga limite (maximo de P(w)) antes do snap-through
%          .dof_apice_v GDL vertical do apice
%          .areas       vetor de areas das duas barras

if nargin < 1 || isempty(a),    a    = 1000;    end
if nargin < 2 || isempty(h),    h    = 100;     end
if nargin < 3 || isempty(E),    E    = 70000;   end
if nargin < 4 || isempty(A),    A    = 100;     end
if nargin < 5 || isempty(dens), dens = 2.7e-6;  end

caso.nome = 'Trelica rasa de 2 barras (von Mises) - validacao analitica';

caso.a = a;  caso.h = h;  caso.A = A;
caso.E = E;  caso.dens = dens;

% -------------------------------------------------------------------------
% GEOMETRIA E TOPOLOGIA
% -------------------------------------------------------------------------
%   no 1 = apoio esquerdo (-a, 0)
%   no 2 = apoio direito  (+a, 0)
%   no 3 = apice          ( 0, h)
caso.nodes0 = [ -a,  a,  0 ;
                 0,  0,  h ];

%   barra 1: no 1 -> no 3
%   barra 2: no 2 -> no 3
caso.elements = [1 2 ;
                 3 3];

caso.areas = [A, A];

% -------------------------------------------------------------------------
% CONDICOES DE CONTORNO
% -------------------------------------------------------------------------
% GDLs: no k ocupa (2k-1) em x e (2k) em y.
%   no 1 -> 1, 2   (fixo)
%   no 2 -> 3, 4   (fixo)
%   no 3 -> 5, 6   (livre)
caso.apoios      = [1 2 3 4];
caso.dof_apice_h = 5;
caso.dof_apice_v = 6;

% -------------------------------------------------------------------------
% CARREGAMENTO (default zero; o ensaio define a carga)
% -------------------------------------------------------------------------
caso.F_total = zeros(6, 1);

% -------------------------------------------------------------------------
% ORACULO ANALITICO
% -------------------------------------------------------------------------
L0 = sqrt(a^2 + h^2);
caso.L0 = L0;

% P(w) conforme a equacao (*) acima. w positivo para BAIXO, P positivo para BAIXO.
caso.P_de_w = @(w) -2 * (E * A * (sqrt(a^2 + (h - w).^2) - L0) / L0) ...
                      .* (h - w) ./ sqrt(a^2 + (h - w).^2);

% Carga critica: maximo de P(w) no intervalo 0 < w < h (ponto limite).
% Obtida por varredura fina seguida de refinamento — suficiente para
% dimensionar ensaios com folga, sem exigir toolbox de otimizacao.
w_varredura = linspace(0, h, 20001);
P_varredura = caso.P_de_w(w_varredura);
[P_max, idx] = max(P_varredura);
caso.P_critica = P_max;
caso.w_critico = w_varredura(idx);

% Inversao numerica: dado P, encontra w no ramo ESTAVEL (antes do ponto
% limite), por bisseccao. Nao usa fzero para funcionar igual em MATLAB e
% Octave sem depender de toolbox.
caso.w_de_P = @(P) inverter_P(P, caso.P_de_w, caso.w_critico);

end


function w = inverter_P(P_alvo, P_de_w, w_critico)
% Bisseccao no ramo estavel [0, w_critico], onde P(w) e crescente.

assert(P_alvo >= 0, 'problema_trelica_rasa_2barras:cargaNegativa', ...
    'Este oraculo cobre apenas cargas nao negativas (para baixo).');

P_max = P_de_w(w_critico);
assert(P_alvo <= P_max, 'problema_trelica_rasa_2barras:cargaAcimaDaCritica', ...
    ['Carga %.6g excede a carga limite %.6g da trelica rasa. ' ...
     'Acima do ponto limite nao existe equilibrio no ramo estavel ' ...
     '(snap-through).'], P_alvo, P_max);

lo = 0;  hi = w_critico;
for k = 1:200
    mid = 0.5 * (lo + hi);
    if P_de_w(mid) < P_alvo
        lo = mid;
    else
        hi = mid;
    end
end
w = 0.5 * (lo + hi);
end
