# Otimização Estrutural com PSO-RID e Análise Não Linear

Implementação do **PSO-RID** (Particle Swarm Optimization Real-Inteiro-Discreto) aplicado
à otimização de peso de treliças com **não linearidade geométrica** ($K_T = K_E + K_G$),
Newton-Raphson incremental e tratamento de restrições pela regra de viabilidade de Deb.

---

## Referências teóricas (fundamentação do código)

Cada fórmula no código cita a referência e o número da equação. As três abaixo são a base:

| Sigla | Referência | Papel |
|---|---|---|
| `[DF2011]` | Datta, D.; Figueira, J.R. *A real-integer-discrete-coded particle swarm optimization for design problems*. Applied Soft Computing 11 (2011) 3625–3633 | **Norte teórico do PSO.** Eq. (1)–(9), Tabelas 1–2 |
| `[DEB2000]` | Deb, K. *An efficient constraint handling method for genetic algorithms*. CMAME 186 (2000) 311–338 | Regra de viabilidade sem parâmetro de penalização (pág. 316) |
| `[HA2003]` | Hadi, M.N.S.; Alvani, K.S. *Discrete Optimum Design of Geometrically Non-Linear Trusses using Genetic Algorithms*. Civil-Comp Press, 2003, Paper 37 | Benchmark de 10 barras; Eq. (13)–(14)–(18) do MEF não linear |

Apoio: Crisfield Vol. 1 (1991) — formulação de base citada por `[HA2003]`.
PDFs em `00_docs/artigos/` e `00_docs/livros/`.

---

## Organização em 3 blocos

O projeto é dividido em três blocos **testáveis isoladamente**:

```text
CPIO III/
│
├── 01_pso_rid/                  ◄── BLOCO 1: otimizador
│   ├── pso_rid.m                    ARQUIVO ÚNICO: loop do PSO + todas as peças
│   │                                como funções locais (ver tabela abaixo)
│   └── pso_classico.m               PSO contínuo padrão (referência)
│
├── 02_fem_nao_linear/           ◄── BLOCO 2: solvers
│   ├── fem_nao_linear_solver.m      Newton-Raphson, K_T = K_E + K_G
│   └── fem_linear_solver.m          análise elástica linear
│
├── 03_orquestrador/             ◄── BLOCO 3: experimentos + definição dos casos
│   ├── main_hadi_nao_linear.m       benchmark principal + DEFINE o caso Hadi
│   ├── main_hadi_linear.m           contraparte linear (comparação)
│   ├── main_awruch_discreto.m       catálogos por barra + DEFINE o caso Awruch
│   ├── main_estudo_estatistico.m    comparação multi-semente de configurações
│   └── auxiliares/                  só blocos longos de impressão/gráfico
│       ├── relatorio_comparativo.m
│       ├── plot_convergencia.m
│       └── salvar_figura.m
│
├── 04_resultados/               figuras e logs gerados
├── 05_legado/                   versões anteriores, para comparação
├── 06_testes/                   matlab.unittest (um arquivo por bloco)
│   └── problema_trelica_rasa_2barras.m   oráculo analítico (fixture de teste)
├── 00_docs/                     artigos, livros, notas e relatórios
└── setup_paths.m
```

### Onde vive cada caso de estudo

Os casos não são arquivos próprios: vivem como função local dentro do orquestrador
que os usa. Como função local não é visível de fora, cada um é exposto por um **acessor**
que devolve só o struct do problema, sem rodar otimização:

| Caso | Definido em | Acessor |
|---|---|---|
| Hadi 10 barras `[HA2003]` | `main_hadi_nao_linear.m` | `main_hadi_nao_linear('caso')` |
| Awruch (catálogos por barra) | `main_awruch_discreto.m` | `main_awruch_discreto('caso')` |
| Treliça rasa (oráculo analítico) | `06_testes/problema_trelica_rasa_2barras.m` | chamada direta |

**O caso Hadi tem fonte única.** Ele é usado por três orquestradores
(`main_hadi_nao_linear`, `main_hadi_linear`, `main_estudo_estatistico`), pelo caso Awruch
— que herda dele geometria, material e cargas — e pelos testes do Bloco 2. Todos passam
pelo mesmo acessor, e não por cópias.

A treliça rasa ficou em `06_testes/` porque **não é um caso de estudo**: é o oráculo com
solução fechada usado só para validar o Newton-Raphson, sem nenhum orquestrador.

### Dentro de `pso_rid.m`

O algoritmo inteiro vive em **um arquivo só**. Cada fórmula do artigo continua sendo uma função separada — só
que como *função local*, no fim do mesmo `.m`. Para navegar, busque pela marca `>>>`:

| Função local | Equação | Papel |
|---|---|---|
| `rid_mapear_dimensoes` | Eq. (5) | dimensões da partícula |
| `rid_decodificar` | Eq. (6) | binário → real/inteiro/discreto |
| `rid_velocidade_binaria` | Eq. (7)-(8), Tab. 1-2 | velocidade em dimensão binária |
| `rid_mutacao_polinomial` | Eq. (9) | mutação polinomial (variável real) |
| `rid_domina_deb` | `[DEB2000]` pág. 316 | regra de viabilidade de Deb |
| `amostrar_particula` | Sec. 4.1 | inicialização aleatória |
| `get_param` | — | default de hiperparâmetro |

Funções locais não são visíveis fora do arquivo. Para exercitar cada peça **isoladamente**
(é o que os testes do Bloco 1 fazem), o solver expõe os handles:

```matlab
aux = pso_rid('auxiliares');
aux.decodificar(bits, mapa, cfg, 'datta')
aux.domina_deb(custo_novo, viol_novo, custo_ref, viol_ref)
```

---

## Como executar

```matlab
setup_paths                  % configura os caminhos

main_hadi_nao_linear         % benchmark principal (não linear)
main_hadi_linear             % contraparte linear
main_awruch_discreto         % catálogos independentes por barra
main_estudo_estatistico(5)   % comparação multi-semente (5 sementes)
```

Todos aceitam argumentos opcionais e devolvem um struct com os resultados:

```matlab
r = main_hadi_nao_linear(123, 10);   % semente 123, 10 execuções
```

---

## Testes

Framework `matlab.unittest`. **Exige MATLAB** — o GNU Octave não implementa esse
framework.

```bash
matlab -batch "cd('06_testes'); run_todos_testes"
```

```matlab
runtests('06_testes')                              % tudo
runtests('06_testes/TestFemNaoLinear.m')           % um bloco só
```

| Arquivo | Bloco | Cobre |
|---|---|---|
| `TestDecodificadorPSO.m` | 1 | Eq. (5), Eq. (6), estratégias de decodificação discreta |
| `TestVelocidadeBinaria.m` | 1 | Eq. (7)–(8), Tabelas 1–2, probabilidade de mutação |
| `TestMutacaoPolinomial.m` | 1 | Eq. (9) e localidade do índice η |
| `TestRegrasDeb.m` | 1 | os 3 critérios de `[DEB2000]` + propriedades da relação |
| `TestFemLinear.m` | 2 | soluções fechadas, simetria, linearidade, equilíbrio |
| `TestFemNaoLinear.m` | 2 | **validação analítica exata**, equilíbrio, limite linear |
| `TestIntegracao.m` | 1+2+3 | interfaces entre blocos, ciclo completo, reprodutibilidade |

---

## Estado da validação

**Newton-Raphson validado contra solução analítica exata.** A treliça rasa de von Mises
(`problema_trelica_rasa_2barras.m`) tem solução fechada; o solver reproduz o deslocamento
do ápice com erro relativo **~1e-11** de 10% a 90% da carga crítica — regime em que a
análise linear erra 37,9%.

> **Nota metodológica importante:** a fórmula analítica foi derivada para **deformação
> de engenharia**, que é a medida usada pelo solver. Fórmulas clássicas de treliça abatida
> (Crisfield, Yaw) usam **deformação de Green** e *não* são diretamente comparáveis — usá-las
> produziria um desacordo que seria erroneamente lido como bug.

### Discrepância do peso — **resolvida** (2026-08-31)

`[HA2003]` enuncia o problema em SI, mas todos os valores são conversões de números
redondos **imperiais**: 9144 mm = 360 in, 50,80 mm = 2 in, 172,25 MPa = 25 ksi,
6,89×10⁴ MPa = 10⁴ ksi, 2770 kg/m³ = 0,1 lb/in³, 445,4 kN = 100 kip, e o catálogo é
`0.1, 1, 2, 5, 8, 12, 15, 18, 20, 25, 30, 35, 45 in²`.

A Tabela 1 foi calculada com **ρ = 0,1 lb/in³ exato (2767,990 kg/m³)** — não com o
2770 kg/m³ que o próprio texto declara:

| Densidade usada | Peso da solução de referência | vs. publicado (2325,2 kg) |
|---|---|---|
| 2770 kg/m³ (texto do artigo) | 2326,942 kg | +0,075% |
| 2767,990 kg/m³ (0,1 lb/in³) | **2325,254 kg** | **+0,002%** |

O resíduo de 0,002% é o arredondamento de impressão do próprio artigo. A explicação vale
para **todas** as colunas da Tabela 1 que listam áreas — inclusive as de
Gutkowski-Zawidzka, obtidas por outro método: as cinco reproduzem o peso publicado com
erro ≤ 0,002%. Travado por `TestFemNaoLinear/testHadi_TodasAsSolucoesPublicadasReproduzemOPeso`.

Os valores exatos estão em `caso.dens_exata` e `caso.E_exato`. **O padrão continua sendo
os valores declarados no texto**, para não invalidar os resultados já registrados em
`04_resultados/`.

### Discrepância do deslocamento — **caracterizada**, resíduo de 0,36%

A solução de referência do artigo excede o limite de deslocamento nesta reimplementação.
Usar o E exato (10⁴ ksi) reduz o excesso de +0,43% para **+0,36%**. O resíduo **não** é:

| Hipótese testada | Resultado |
|---|---|
| Erro de convergência do Newton-Raphson | Descartada — `max\|u\|` idêntico até a 6ª casa para `n_inc` de 1 a 200 |
| Numeração das diagonais (barras 7–10, mesmo comprimento) | Descartada — a ordem adotada dá o **menor** deslocamento das 24 permutações; nenhuma torna a referência viável |
| Arredondamento das áreas (mm² vs in² exato) | Descartada — diferença de 0,0008 mm |
| Interpretação da restrição | Descartada — `[HA2003]` Eq. (2) restringe **todos** os GDLs, igual ao código |

**O que sobra é a natureza do problema:** ele é governado por deslocamento, e o ótimo fica
*sobre* a restrição. As seis soluções publicadas ficam a menos de 1% do limite (margens de
−0,93% a +0,36%); uma diferença de formulação de 0,36% basta para virar a viabilidade.
O solver reproduz corretamente o estado de restrição das outras cinco — todas apenas
viáveis. Sob este modelo, a solução publicada mais leve que é *viável* é a G-Z Enumeration
do Case 2, com 2339,9 kg a 50,708 mm.

---

## Decisões de projeto (D1–D8)

O solver diverge de `[DF2011]`/`[HA2003]` em pontos catalogados como `D1`–`D8`. Cada um é uma decisão
explícita, com justificativa no código e evidência no relatório. Detalhes completos em
`00_docs/notas_e_relatorios/relatorio_decisoes_2026-08-31.md`.

| | Tema | Decisão | Executa? |
|---|---|---|---|
| **D1** | Decodificação discreta | **`'proporcional'`** (desvia do artigo) | sim, sempre |
| D2 | Decodificação inteira | fiel ao artigo | não — todos os casos são tipo `'D'` |
| D3 | Velocidade binária assimétrica | fiel ao artigo | sim, sempre |
| D4 | Mutação binária `pm` | fiel ao artigo | sim, sempre |
| D5 | Mutação polinomial (variável real) | fiel ao artigo | não — todos os casos são tipo `'D'` |
| D6 | Reinicialização de partículas | **em aberto** — heurística fora do artigo | sim, sempre |
| D7 | Escala da violação de estouro | violação = `Inf` | só sob `'datta'` |
| **D8** | Medida de violação | **normalizada linear, `[HA2003]` Eq. (9)** | sim, sempre |

### D1 — por que a decodificação não segue o artigo

`[DF2011]` prescreve tratar índices fora do catálogo como restrição (Eq. 12 e Sec. 5.2:
*"its integer limits are set automatically as [1,42]"*). O problema é a **escala**: com
N=13 opções em B=4 bits sobram 3 dos 16 códigos por variável, e com 10 variáveis isso
compõe para `(13/16)¹⁰ = 12,5%` — só uma partícula em oito é inteiramente válida. O artigo
nunca testou esse regime; no seu único exemplo com discreta há **uma** variável.

Evidência (30 sementes pareadas, orçamento igualado em avaliações FEM, Hadi linear):

| Configuração | Melhor | Média | Desvio |
|---|---|---|---|
| `'datta'` | 2357,97 | 2572,64 | 173,05 |
| **`'proporcional'`** | **2342,11** | **2504,72** | **97,40** |

Variância: **F = 3,16 (crítico 2,10) → significativo.** Média: t = −1,836 (sem evidência),
mas o teste do sinal dá 21/30 (**p = 0,043, significativo**) — o t perde pela cauda pesada
das diferenças.

Custo assumido: viés estatístico. As áreas 65, 5161 e 12903 mm² recebem 2 dos 16 códigos
cada, contra 1 das outras dez. O modo `'datta'` segue disponível por parâmetro.

Sob `'datta'`, os códigos que caem fora do catálogo recebem **violação = `Inf`**: uma
partícula fora do catálogo não é um projeto ruim, é a ausência de projeto, e deve perder
para qualquer projeto avaliável. `Inf` e não uma constante grande porque violações reais
neste problema chegam a 4,5×10⁸ — qualquer constante fixa seria ultrapassada.

### D8 — por que a violação é normalizada

`avaliar_projeto` somava os quadrados dos excessos brutos, ou seja **MPa² com mm²**. A razão
entre as duas parcelas passava a depender da unidade escolhida, não da física — e o critério
3 de `[DEB2000]` ordena inviáveis **somente** pela violação. `[HA2003]` Eq. (9) já prescreve
a forma adimensional, que o código agora usa:

    g_σ,i = |σ_i|/σ_adm − 1 ≤ 0        g_d,j = |d_j|/d_adm − 1 ≤ 0
    violação = Σ max(0, g)

Medido em 2000 pares de projetos inviáveis: a medida antiga **inverte a ordenação em 1,7%
dos pares** quando os deslocamentos passam de mm para m; a nova é invariante por construção.

O efeito no desempenho é **nulo** — 30 sementes pareadas, mesmo ótimo (2342,11), mesma
viabilidade (30/30), t = −0,857 e F = 0,918, ambos sem evidência. A mudança é de correção,
não de desempenho, e custa zero.

### Mecanismos fora da bibliografia que foram avaliados e mantidos

| Mecanismo | Veredito |
|---|---|
| `V_max = (x_max−x_min)/2` nas dimensões reais | mantido — prática padrão do PSO desde Shi & Eberhart (1998) |
| Clamp da posição real em `[min, max]` | mantido — **obrigatório**: sem ele `mutacao_polinomial` retorna complexo silenciosamente |
| Limites de índice / `viol_estrutural` | mantido — inerte sob `'proporcional'`, mas sustenta o tipo `'I'` e o modo `'datta'` |

Os dois primeiros só executam com variáveis do tipo `'R'`, que nenhum estudo de caso atual
usa. O terceiro fica em zero sob o padrão, mas removê-lo custaria o tipo inteiro da Eq. (12)
de `[DF2011]` e a reprodutibilidade da evidência de D1.

### Q2 — respondida, hipótese refutada

A suspeita de que a auto-adaptação de hiperparâmetros de `[DF2011]` prejudicava a
convergência era artefato de amostra de tamanho 1: com 10 sementes pareadas, +31,30 kg,
t = 0,631, 4/10 vitórias. **Sem evidência.**

**Ainda em aberto:** o estudo pareado só foi feito com FEM **linear** (~5 min); a
contraparte não linear custa ~4 h e não foi executada, então as conclusões de D1 podem não
transferir. E D6 (reinicialização) nunca foi medido.

---

## Documentação

`00_docs/notas_e_relatorios/`:

- `plano_final_reorganizacao_2026-08-28.md` — plano dos 3 blocos e referências
- `relatorio_mudancas_2026-08-28.md` — o que mudou, por quê, com evidência numérica
- `backup_contexto_2026-08-28.md` — análise de validação teórica original
- **`relatorio_decisoes_2026-08-31.md`** — decisão D1, resolução da discrepância de peso
  (unidades imperiais), caracterização da discrepância de deslocamento, e as pendências
  menores fechadas
