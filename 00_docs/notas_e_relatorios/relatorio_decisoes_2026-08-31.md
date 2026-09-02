# Relatório de decisões — 2026-08-31

> Continuação de `relatorio_mudancas_2026-08-28.md`. Este documento fecha as pendências
> que aquele deixou em aberto e registra as decisões **D1** e **D7** com a evidência que
> as sustenta. Tudo aqui é reprodutível: os testes citados travam cada resultado, e os
> logs estão em `04_resultados/logs/`.

---

## Sumário das decisões

| Item | Estado antes | Decisão | Base |
|---|---|---|---|
| **D1** decodificação discreta | `'datta'` padrão, `'proporcional'` disponível | **`'proporcional'` vira o padrão** | 30 sementes pareadas, teste F significativo |
| **D7** escala da violação de estouro | não catalogado (defeito latente) | **`'dominante'` vira o padrão** | correção semântica; efeito nulo medido |
| Discrepância do peso vs `[HA2003]` | documentada, inexplicada | **resolvida** | problema é imperial; ρ = 0,1 lb/in³ |
| Discrepância do deslocamento | documentada, inexplicada | **caracterizada**, resíduo 0,36% | 4 hipóteses testadas e descartadas |
| Figuras linear ≡ não linear | suspeita de bug | **falso alarme** | valores extraídos dos `.fig` diferem |
| Monotonicidade do `cost_history` | testada como não crescente | **premissa corrigida** | Deb crit. 1 permite o custo subir |
| **D8** medida de violação | quadrática, sem normalizar | **normalizada linear, `[HA2003]` Eq. (9)** | correção de unidades; efeito nulo medido |
| Limite de velocidade real (`V_max`) | não catalogado | **mantido**, agora com citação | Shi & Eberhart (1998) |
| Clamp da posição real | não catalogado | **mantido**, é obrigatório | sem ele `mutacao_polinomial` retorna complexo |
| Limites de índice / `viol_estrutural` | suspeita de código morto | **mantido**, dormente sob o padrão | sustenta tipo `'I'` e modo `'datta'` |

---

## 1. A discrepância do peso — resolvida

`[HA2003]` enuncia o problema em SI, mas **todos** os seus valores são conversões de
números redondos imperiais:

| SI declarado no artigo | Imperial exato |
|---|---|
| 9144 mm | 360 in |
| 50,80 mm | 2 in |
| 172,25 MPa | 25 ksi |
| 6,89×10⁴ MPa | 10⁴ ksi = 68947,573 MPa |
| 2770 kg/m³ | 0,1 lb/in³ = 2767,990 kg/m³ |
| 445,4 kN | 100 kip = 444,822 kN |
| catálogo em mm² | 0.1, 1, 2, 5, 8, 12, 15, 18, 20, 25, 30, 35, 45 in² |

A Tabela 1 foi calculada com **ρ = 0,1 lb/in³ exato**, não com o 2770 kg/m³ que o próprio
texto declara:

| Densidade | Peso da solução de referência | vs. 2325,2 kg publicado |
|---|---|---|
| 2770 kg/m³ (texto) | 2326,942 kg | +0,075% |
| 2767,990 kg/m³ (0,1 lb/in³) | **2325,254 kg** | **+0,002%** |

O resíduo de 0,002% é o arredondamento de impressão do artigo (uma casa decimal).

**A explicação vale para toda a Tabela 1.** Aplicada às cinco colunas que listam áreas —
incluindo as de Gutkowski-Zawidzka, obtidas por outro método — todas reproduzem o peso
publicado com erro ≤ 0,002%:

```
C1 G-Z Enumeration  : 2429,51 kg (publicado 2429,5, +0,0006%)
C1 Hadi discreto    : 2423,21 kg (publicado 2423,2, +0,0002%)
C2 G-Z Sequential   : 2340,41 kg (publicado 2340,4, +0,0004%)
C2 G-Z Enumeration  : 2339,93 kg (publicado 2339,9, +0,0015%)
C2 Hadi discreto    : 2325,25 kg (publicado 2325,2, +0,0023%)
```

**Onde ficou:** `caso.dens_exata` e `caso.E_exato`, na função local `caso_hadi_10barras`
de `main_hadi_nao_linear.m`. **O padrão continua sendo os valores do texto do artigo**,
para não invalidar os resultados já registrados em `04_resultados/`.

**Travado por:** `TestFemNaoLinear/testHadi_UnidadesImperiaisExplicamOPeso` e
`testHadi_TodasAsSolucoesPublicadasReproduzemOPeso`.

---

## 2. A discrepância do deslocamento — caracterizada

A solução de referência do artigo excede o limite de deslocamento nesta reimplementação.
Usar o E exato reduz o excesso de +0,43% para **+0,36%**. O resíduo foi investigado
hipótese a hipótese:

| Hipótese | Teste | Resultado |
|---|---|---|
| Erro de convergência do Newton-Raphson | `n_inc` de 1 a 200 | `max\|u\|` idêntico até a 6ª casa — **descartada** |
| Numeração das diagonais (barras 7–10 têm o mesmo comprimento, então permutá-las não muda o peso — só o deslocamento) | todas as 24 permutações | a ordem adotada dá o **menor** deslocamento (50,984); a seguinte é 51,025 e as demais passam de 68 mm; nenhuma torna a referência viável — **descartada** |
| Arredondamento das áreas (mm² publicado vs in² exato) | ambas | diferença de 0,0008 mm — **descartada** |
| Interpretação da restrição | `[HA2003]` Eq. (2) | `d_j ≤ d_j^a, j = 1..m` restringe **todos** os GDLs, igual ao código — **descartada** |

**O que sobra é a natureza do problema.** Ele é governado por deslocamento e o ótimo fica
*sobre* a restrição. As seis soluções publicadas ficam a menos de 1% do limite:

```
C1 G-Z Seq  :  -6,976%      C2 G-Z Seq  :  -0,636%
C1 G-Z Enum :  -0,509%      C2 G-Z Enum :  -0,182%
C1 Hadi     :  -0,931%      C2 Hadi     :  +0,363%  <-- único acima
```

Uma diferença de formulação de 0,36% basta para virar a viabilidade. O solver reproduz
corretamente o estado de restrição das outras cinco (todas apenas viáveis), o que é a
evidência de que o modelo está certo e a discrepância é específica daquele projeto.

Sob este modelo, a solução publicada mais leve que é **viável** é a G-Z Enumeration do
Case 2: 2339,9 kg a 50,708 mm.

**Travado por:** `TestFemNaoLinear/testHadi_NumeracaoDasDiagonaisEhAMelhorPossivel`.

---

## 3. [D1] Decodificação discreta — `'proporcional'` vira o padrão

### O que o artigo diz sobre códigos fora de faixa

`[DF2011]` é silencioso na formulação, mas responde pelos exemplos numéricos. No trem de
engrenagens (Eq. 12) escreve os limites do inteiro como **restrição explícita**:

> `subject to: g(x) ≡ 12 ≤ za, zb, zc, zd ≤ 60`

Com 6 bits, `z` varre 0–63; os códigos 0–11 e 61–63 viram violação de restrição, tratada
pela regra de Deb. E na mola (Sec. 5.2) estende isso ao discreto:

> "Since d is a discrete variable with 42 allowable values given in Table 5, **its integer
> limits are set automatically as [1,42]**."

Ou seja: código fora de faixa **não** é descartado, remapeado nem clampado — é uma solução
inviável. Era exatamente o que o modo `'datta'` fazia.

### Por que ainda assim o padrão mudou

**A escala do desperdício é outra.** Na mola, Datta tem *uma* variável discreta (42 em 64 =
34% de códigos mortos). O benchmark de Hadi tem **dez**, e o desperdício compõe:

```
(13/16)^10 = 12,5%
```

Apenas uma partícula em oito é inteiramente válida; **87,5% nunca viram projeto avaliável**
(medido: 86,9% em 4000 amostras). `[DF2011]` nunca testou esse regime.

### Evidência

30 sementes pareadas, orçamento igualado em 15.000 avaliações FEM, benchmark Hadi linear.
Log: `04_resultados/logs/estudo_estatistico_linear_20260831_143039.txt`

| Configuração | Melhor | Média | Desvio |
|---|---|---|---|
| `'datta'` + violação literal | 2357,97 | 2572,64 | 173,05 |
| `'datta'` + dominante [D7] | 2342,11 | 2579,14 | 214,66 |
| **`'proporcional'`** | **2342,11** | **2504,72** | **97,40** |

| Teste | Estatística | Veredito |
|---|---|---|
| F nas variâncias, `proporcional` vs `datta` | F = 3,16 (crítico 2,10) | **significativo** |
| F nas variâncias, `proporcional` vs `datta+dominante` | F = 4,86 (crítico 2,10) | **significativo** |
| t pareado nas médias, `proporcional` − `datta` | t = −1,836 (crítico 2,045) | sem evidência |
| Teste do sinal, mesma comparação | 21/30, p = 0,043 | **significativo** |

O t pareado perde por causa da cauda pesada das diferenças — ressalva que o próprio estudo
já registrava. O teste do sinal, livre de distribuição, pega.

**Mecanismo, não sorte:** sem códigos mortos, 100% das partículas viram projeto avaliável
contra 12,5% no `'datta'`. O enxame carrega cerca de oito vezes mais informação de projeto
por iteração.

### Custo assumido

O `'proporcional'` tem viés estatístico. Com N=13 em B=4 bits, os 16 códigos se distribuem
assim:

```
índice:  1  2  3  4  5  6  7  8  9 10 11 12 13
códigos: 2  1  1  1  2  1  1  1  2  1  1  1  1
```

As áreas 65, 5161 e 12903 mm² recebem o dobro de chance das outras dez. É um desvio
declarado de `[DF2011]`, assumido em troca da consistência medida acima. O modo `'datta'`
continua disponível por parâmetro para reproduzir o comportamento do artigo.

---

## 4. [D7] Escala da violação de estouro — `'dominante'` vira o padrão

### O defeito

Quando uma partícula estoura o catálogo, o código atribuía `violação = excesso de índice`
(1, 2, 3…). Mas a violação de um projeto real vem em MPa² + mm². As escalas não são
comensuráveis, e a regra de Deb compara as duas diretamente no critério 3:

| | Faixa | Mediana |
|---|---|---|
| Violação estrutural (índices) | 1 – 15 | 4 |
| Violação FEM (MPa² + mm²) | 2×10⁻⁵ – 4,5×10⁸ | 9,1×10⁴ |

Medido em 4000 amostras: **em 97,9% dos pares, Deb prefere a partícula fora do catálogo a
um projeto real.** O efeito é uma armadilha de memória — um `pbest` que é um estouro com
violação 1 só aceita ser trocado por uma solução plenamente viável ou por outra inviável
com violação < 1, o que é impossível (a violação estrutural é inteira ≥ 1). Ele congela
numa posição sem significado físico.

### A correção

Uma partícula fora do catálogo não é um projeto ruim: é a **ausência** de projeto. Deve
perder para qualquer projeto avaliável. Isso é uma quarta regra somada às três de Deb:

```
1. viável                     vence  qualquer inviável         [DEB2000]
2. entre viáveis              vence  menor custo               [DEB2000]
3. entre inviáveis avaliáveis vence  menor violação            [DEB2000]
4. avaliável                  vence  estruturalmente inviável  [D7]
```

Implementada como `violação = Inf`, o que faz as três regras originais produzirem a quarta
sozinhas: `Inf < Inf` é falso (nenhuma preferência entre estouros) e `1e8 < Inf` é
verdadeiro (avaliável vence estouro).

**Sem número mágico, deliberadamente.** Uma constante grande falharia em silêncio: foram
medidas violações reais de até 4,5×10⁸. A comparação tem de ser por classe, não por
magnitude.

**Sem ordenação entre estouros**, também deliberadamente. A pressão que traz a partícula de
volta ao catálogo vem do g-best pela Eq. (1), não da comparação entre posições sem
significado físico.

### Por que virou padrão mesmo sem ganho

O efeito no resultado é **nulo**: +6,51 kg, t = 0,169, 17/30 vitórias — ruído puro.

A troca é de **correção semântica**, não de desempenho. E não custa fidelidade: `[DF2011]`
prescreve tratar os limites do índice como restrição, mas é **silencioso sobre a escala da
violação**. Usar o excesso bruto foi invenção desta implementação, não leitura do artigo.

> **Nota metodológica.** A hipótese original era que corrigir isso reduziria a dispersão do
> `'datta'` (173,5 vs 88,8 do proporcional). **Foi refutada:** a dispersão subiu para 214,7.
> Explicação plausível: com `'dominante'`, `pbest` e `gbest` só podem ser projetos
> avaliáveis, então os 87% de partículas em estouro ficam com `pbest` congelado na posição
> inicial aleatória — antes ao menos rastreavam o estouro menos ruim. Removeu-se um atrator
> espúrio e, junto, o único sinal de `pbest` que a maioria do enxame tinha.

**Travado por:** `TestRegrasDeb/testD7_*` (quatro testes, incluindo um de regressão que
documenta a patologia do modo literal).

---

## 5. Pendências menores fechadas

### Figuras linear ≡ não linear — falso alarme

`backup_contexto_2026-08-28.md` item 5 suspeitava de bug na separação dos casos, porque as
duas figuras pareciam mostrar 2342,11 kg. Extraindo os dados dos `.fig` (via `openfig` +
`YData`, em vez da leitura visual do `.png`, que foi a origem do engano):

| Figura | Melhor valor plotado |
|---|---|
| `convergencia_hadi_linear.fig` | 2342,1077 kg |
| `convergencia_hadi_nao_linear.fig` | 2340,4613 kg |

São **diferentes**. Não há problema na separação dos casos.

### `cost_history` não é monotônico — premissa do teste corrigida

`TestIntegracao` exigia que o histórico de custo do g-best fosse não crescente. **Isso está
errado sob Deb:** pelo critério 1, um viável é preferido a qualquer inviável
independentemente do custo, então o custo **sobe** quando o enxame troca um inviável-barato
por um viável-caro. Observado com semente 42: 3916,83 → 4543,62 na iteração 1→2.

O teste passava por acidente. Com `'datta'`, as primeiras iterações eram dominadas por
estouros com custo `Inf`, que o filtro `isfinite` removia — escondendo a fase inviável.
Trocar o padrão para `'proporcional'` expôs a premissa errada.

O invariante de fato garantido pelos critérios 1 e 3 é sobre a **violação**, que nunca
piora. Foi adicionado `details.viol_history` ao `pso_rid`, e o teste agora verifica:

- a violação do g-best é monotonicamente não crescente;
- atingida a viabilidade, ela nunca se perde;
- **dentro** da região viável, aí sim o custo só pode cair (critério 2).

---

## 6. [D8] Medida de violação — normalizada, `[HA2003]` Eq. (9)

### O defeito

`avaliar_projeto` somava os quadrados dos excessos **brutos**:

```matlab
violacao = sum((|sigma| - sigma_max).^2) + sum((|u| - d_max).^2);
```

Isso soma **MPa² com mm²**. O número resultante não tem unidade coerente, e a razão entre
as duas parcelas passa a depender de em que unidade o problema foi escrito, não da física.
Com σ_adm = 172,25 MPa e d_adm = 50,80 mm, 1 mm de excesso de deslocamento pesava o mesmo
que 1 MPa de excesso de tensão — embora o primeiro consuma 2,0% da folga e o segundo 0,6%.

Isso importa porque o **critério 3 de `[DEB2000]` ordena inviáveis SOMENTE pela violação**.
Uma medida incomensurável faz a busca ordenar pela unidade escolhida.

`[HA2003]` Eq. (9) já prescreve a forma correta, e o código não a seguia:

    g_σ,i = |σ_i| / σ_adm − 1 ≤ 0
    g_d,j = |d_j| / d_adm − 1 ≤ 0

### A correção

```matlab
g_sigma = abs(Sigma(:))   / caso.sigma_max - 1;
g_desl  = abs(u_livre(:)) / caso.d_max     - 1;
violacao = sum(max(0, g_sigma)) + sum(max(0, g_desl));
```

Linear, não quadrática: o quadrado **amplifica** a incomensurabilidade (MPa² contra mm²) em
vez de corrigi-la, e não consta nem de `[HA2003]` nem de `[DEB2000]`.

Aplicada nas **cinco** cópias da medida — `avaliar_projeto` dos quatro orquestradores mais
`custo_linear` de `TestIntegracao.m`, que antes divergia silenciosamente do que roda de fato.

### Quanto o defeito de fato mordia

Amostragem de 4000 projetos inviáveis do catálogo, comparados dois a dois (2000 pares):

| Medida | Ordenação inverte ao trocar mm → m | Faixa observada |
|---|---|---|
| quadrática bruta | **34 de 2000 pares (1,7%)** | 4,07×10⁻⁴ a 5,24×10⁸ |
| normalizada | **0 — invariante por construção** | 3,97×10⁻⁴ a 7,00×10² |

As duas discordam da ordenação em 3,6% dos pares. O número é baixo porque neste problema
o deslocamento domina: a parcela de tensão raramente arbitra. **Num problema em que as duas
famílias competissem de verdade, o efeito seria muito maior.**

### Efeito no desempenho: nulo

Comparação pareada, 30 sementes, orçamento de 15000 avaliações, Hadi linear:

| Medida | Melhor | Média | Desvio | Pior | Viáveis |
|---|---|---|---|---|---|
| quadrática bruta | 2342,11 | 2502,77 | 97,81 | 2760,90 | 30/30 |
| **normalizada** | **2342,11** | 2522,01 | 102,10 | **2699,46** | 30/30 |

| Teste | Estatística | Veredito |
|---|---|---|
| t pareado nas médias | t = −0,857 (crítico 2,045) | sem evidência |
| Teste do sinal | 10/30, p = 0,099 | sem evidência |
| F nas variâncias | F = 0,918 (crítico 2,101) | sem evidência |

**A mudança não melhora o resultado — e não é por isso que foi feita.** Ambas atingem o
mesmo ótimo (2342,11) e a mesma viabilidade (30/30). A justificativa é de correção: a
medida anterior era dependente de unidade e divergia de `[HA2003]`. O custo medido de
adotar a forma correta é zero.

---

## 7. `V_max`, clamp da posição real e limites de índice — os três mantidos

Três mecanismos não descritos em `[DF2011]` foram avaliados para remoção. Todos ficam,
por motivos diferentes.

### Limite de velocidade real — mantido, agora citado

`v_lim = (x_max − x_min)/2`. Não está em `[DF2011]`, mas **é prática padrão do PSO contínuo**
desde Shi & Eberhart (1998), *A modified particle swarm optimizer*, IEEE ICEC, p. 69–73:
sem teto, o termo de inércia `w·v` acumula e a partícula salta de um extremo do domínio ao
outro, perdendo a capacidade de refinar. A recomendação usual é `V_max = k·(x_max − x_min)`
com k ∈ [0,1; 1,0]; k = 0,5 está no meio da faixa. O comentário no código foi trocado de
"estabilização numérica" (vago) para a citação.

### Clamp da posição real — mantido, é obrigatório

Não é opcional. Sem ele, uma dimensão real fora de `[min, max]` quebra dois consumidores:

1. **`mutacao_polinomial` retorna número complexo.** Ela eleva `(x_u − x)/(x_u − x_l)` à
   potência `(η+1)` com η **real**. Base negativa com expoente fracionário dá resultado
   complexo — e o `max/min` final daquela função **não detecta**, porque MATLAB compara
   complexos por módulo. A corrupção seria silenciosa.
2. **`decodificar` repassa o valor cru à função de custo.** Uma área negativa produz matriz
   de rigidez sem sentido físico.

`[DF2011]` trata limites como restrição no caso **inteiro** (ver D2), mas ali o valor cru
ainda é avaliável. Aqui não é, então a saturação é a única saída coerente.

### Limites de índice / `viol_estrutural` — mantidos, dormentes

A suspeita era que a decodificação `'proporcional'` tivesse tornado a máquina de estouro
código morto. Ela a tornou **inerte nos casos atuais, não morta**. Há três caminhos:

| Caminho | Ainda gera `viol_estrutural`? |
|---|---|
| discreto `'proporcional'` (padrão) | **não** — `idx = min(max(idx,1),N)` satura por construção |
| discreto `'datta'` | **sim** — índice > N |
| inteiro `'I'` | **sim, sempre** — fora de `[min, max]` |

Como os três estudos de caso usam só variáveis `'D'` com o padrão, na prática
`viol_estrutural` e `n_overflow` ficam sempre em zero neles. Remover o mecanismo custaria:

- o **tipo inteiro**, exigido pela Eq. (12) de `[DF2011]` (trem de engrenagens, 12 ≤ z ≤ 60
  em 6 bits → 0..63) e exercitado em 6 pontos de `TestDecodificadorPSO`;
- o **modo `'datta'`**, usado como configuração de comparação em `main_estudo_estatistico.m`
  — a evidência que decidiu D1 deixaria de ser reproduzível.

Foi adicionado ao cabeçalho de `decodificar` um bloco "QUANDO ISSO OCORRE", para que o
mecanismo não seja procurado num caminho que nunca executa.

---

## Estado do catálogo D1–D7

| | Tema | Executa nos casos atuais? | Estado |
|---|---|---|---|
| D1 | Decodificação discreta | sim, sempre | **decidido**: `'proporcional'` |
| D2 | Decodificação inteira | **não** (todos os casos são tipo `'D'`) | fiel ao artigo, sem uso prático |
| D3 | Velocidade binária assimétrica | sim, sempre | fiel ao artigo, travado por teste |
| D4 | Mutação binária `pm` | sim, sempre | fiel ao artigo, travado por teste |
| D5 | Mutação polinomial (variável real) | **não** (todos os casos são tipo `'D'`) | fiel ao artigo, sem uso prático |
| D7 | Escala da violação de estouro | só sob `'datta'` | **decidido**: `'dominante'` |
| D8 | Medida de violação | sim, sempre | **decidido**: normalizada, `[HA2003]` Eq. (9) |

**D6 ainda em aberto:** a reinicialização de partículas estagnadas
(`reinit_freq` / `reinit_pct`) **não existe em `[DF2011]`** — é heurística herdada do
legado, ligada por padrão, e altera a busca em todos os experimentos. Nunca foi catalogada
nem medida. Use `reinit_freq = 0` para o comportamento estrito do artigo.

**Também em aberto:** o estudo pareado só foi feito com FEM **linear**. A contraparte não
linear custa cerca de 4 h e não foi executada; as conclusões de D1 podem não transferir.
