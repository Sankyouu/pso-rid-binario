# CPIO III — Otimização Estrutural com PSO-RID e Análise Não Linear

Implementação do **PSO-RID** (Particle Swarm Optimization Real-Inteiro-Discreto) aplicado
à otimização de peso de treliças com **não linearidade geométrica** ($K_T = K_E + K_G$),
Newton-Raphson incremental e tratamento de restrições pela regra de viabilidade de Deb.

---

## 📚 Referências teóricas (fundamentação do código)

Cada fórmula no código cita a referência e o número da equação. As três abaixo são a base:

| Sigla | Referência | Papel |
|---|---|---|
| `[DF2011]` | Datta, D.; Figueira, J.R. *A real-integer-discrete-coded particle swarm optimization for design problems*. Applied Soft Computing 11 (2011) 3625–3633 | **Norte teórico do PSO.** Eq. (1)–(9), Tabelas 1–2 |
| `[DEB2000]` | Deb, K. *An efficient constraint handling method for genetic algorithms*. CMAME 186 (2000) 311–338 | Regra de viabilidade sem parâmetro de penalização (pág. 316) |
| `[HA2003]` | Hadi, M.N.S.; Alvani, K.S. *Discrete Optimum Design of Geometrically Non-Linear Trusses using Genetic Algorithms*. Civil-Comp Press, 2003, Paper 37 | Benchmark de 10 barras; Eq. (13)–(14)–(18) do MEF não linear |

Apoio: Crisfield Vol. 1 (1991) — formulação de base citada por `[HA2003]`.
PDFs em `00_docs/artigos/` e `00_docs/livros/`.

---

## 🧱 Organização em 3 blocos

O projeto é dividido em três blocos **testáveis isoladamente**:

```text
CPIO III/
│
├── 01_pso_rid/                  ◄── BLOCO 1: otimizador
│   ├── pso_rid.m                    loop do PSO (chama os auxiliares abaixo)
│   ├── rid_mapear_dimensoes.m       Eq. (5)      dimensões da partícula
│   ├── rid_decodificar.m            Eq. (6)      binário → real/inteiro/discreto
│   ├── rid_velocidade_binaria.m     Eq. (7)-(8)  Tabelas 1 e 2
│   ├── rid_mutacao_polinomial.m     Eq. (9)      mutação polinomial
│   ├── rid_domina_deb.m             pág. 316     regra de viabilidade de Deb
│   └── pso_classico.m               PSO contínuo padrão (referência)
│
├── 02_fem_nao_linear/           ◄── BLOCO 2: análise estrutural
│   ├── solver/                      2a — genérico, NÃO conhece nenhum problema
│   │   ├── fem_nao_linear_solver.m      Newton-Raphson, K_T = K_E + K_G
│   │   └── fem_linear_solver.m          análise elástica linear
│   └── problemas/                   2b — SÓ parâmetros, um arquivo por caso
│       ├── problema_hadi_10barras.m          benchmark [HA2003]
│       ├── problema_awruch_10barras.m        catálogos independentes por barra
│       └── problema_trelica_rasa_2barras.m   validação analítica exata
│
├── 03_orquestrador/             ◄── BLOCO 3: experimentos
│   ├── main_hadi_nao_linear.m       benchmark principal
│   ├── main_hadi_linear.m           contraparte linear (comparação)
│   ├── main_awruch_discreto.m       catálogos por barra
│   ├── main_estudo_estatistico.m    comparação multi-semente de configurações
│   └── auxiliares/                  só blocos longos de impressão/gráfico
│       ├── relatorio_comparativo.m
│       ├── plot_convergencia.m
│       └── salvar_figura.m
│
├── 04_resultados/               figuras e logs gerados
├── 05_legado/                   versões anteriores, para comparação
├── 06_testes/                   matlab.unittest (um arquivo por bloco)
├── 00_docs/                     artigos, livros, notas e relatórios
└── setup_paths.m
```

**Por que essa divisão:** o solver do Bloco 2a nunca contém dados de um problema
específico. Adicionar um novo caso de estudo é escrever um arquivo pequeno em
`problemas/` com nós, barras, apoios, cargas e material — sem duplicar solver. O Bloco 3
concentra cada experimento em um único `.m` autocontido, deixando fora apenas os blocos
longos de saída.

---

## 🚀 Como executar

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

## 🧪 Testes

Framework `matlab.unittest`. **Exige MATLAB** — o GNU Octave não implementa esse
framework (o Octave continua servindo para rodar os solvers diretamente).

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

## ✅ Estado da validação

**Newton-Raphson validado contra solução analítica exata.** A treliça rasa de von Mises
(`problema_trelica_rasa_2barras.m`) tem solução fechada; o solver reproduz o deslocamento
do ápice com erro relativo **~1e-11** de 10% a 90% da carga crítica — regime em que a
análise linear erra 37,9%.

> ⚠️ **Nota metodológica importante:** a fórmula analítica foi derivada para **deformação
> de engenharia**, que é a medida usada pelo solver. Fórmulas clássicas de treliça abatida
> (Crisfield, Yaw) usam **deformação de Green** e *não* são diretamente comparáveis — usá-las
> produziria um desacordo que seria erroneamente lido como bug.

**Discrepâncias conhecidas e documentadas** (travadas por testes de caracterização, para
que uma mudança de tamanho falhe o teste em vez de passar despercebida):

| Item | Calculado | Publicado `[HA2003]` | Diferença |
|---|---|---|---|
| Peso da solução de referência | 2326,94 kg | 2325,2 kg | +0,075% |
| Deslocamento máximo | 51,02 mm | limite 50,80 mm | +0,43% |

**Questões em aberto** — ver `00_docs/notas_e_relatorios/relatorio_mudancas_2026-08-28.md`
e responder com `main_estudo_estatistico`:

- Q1: a decodificação discreta literal de `[DF2011]` é melhor que a heurística proporcional?
- Q2: a auto-adaptação de hiperparâmetros prejudica a convergência frente a valores fixos?

Ambas ficaram **inconclusivas** com execuções únicas (diferenças dentro do ruído).

---

## 📄 Documentação

`00_docs/notas_e_relatorios/`:

- `plano_final_reorganizacao_2026-08-28.md` — plano dos 3 blocos e referências
- `relatorio_mudancas_2026-08-28.md` — o que mudou, por quê, com evidência numérica
- `backup_contexto_2026-08-28.md` — análise de validação teórica original
