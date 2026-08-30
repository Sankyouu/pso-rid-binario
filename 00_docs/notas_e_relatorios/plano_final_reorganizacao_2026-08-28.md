# Plano Final — Reorganização do Projeto CPIO III em 3 Blocos

> Documento de fechamento de plano (2026-08-28), reunindo (A) o diagnóstico e as
> referências teóricas levantadas na validação e (B) a proposta de reorganização
> estrutural do projeto em 3 blocos testáveis isoladamente. Este arquivo é a fonte de
> verdade do plano; ver também `backup_contexto_2026-08-28.md` para a cópia de segurança
> apenas da Parte A, e `relatorio_mudancas_2026-08-28.md` para o rastreamento das
> mudanças conforme forem aplicadas.

---

# Parte A — Diagnóstico e referências teóricas

## A.1 Referências-norte confirmadas

| Bloco | Referência primária | Papel |
|---|---|---|
| PSO-RID | Datta, D.; Figueira, J.R. — *A real-integer-discrete-coded particle swarm optimization for design problems*, Applied Soft Computing 11 (2011) 3625–3633 | Fórmulas (Eq. 1–9) e Tabelas 1–2 do encoding binário — **norte teórico do bloco PSO** |
| PSO-RID (restrições) | Deb, K. — *An efficient constraint handling method for genetic algorithms*, Comput. Methods Appl. Mech. Engrg. 186 (2000) 311–338 | Regra de dominância de 3 critérios (pág. 316), já confirmada batendo 1:1 com o código |
| FEM Não-Linear | Hadi, M.N.S.; Alvani, K.S. — *Discrete Optimum Design of Geometrically Non-Linear Trusses using Genetic Algorithms*, AICC2003 | Geometria/material/catálogo/carga do benchmark de 10 barras; Eq. (13)-(14) para K_T = K_E + K_G |
| FEM Não-Linear (teoria de base) | Crisfield, M.A. — *Non-Linear Finite Element Analysis of Solids and Structures*, Vol. 1 (1991) | Formulação corrotacional/incremental citada pelo Hadi [16] — já em `00_docs/livros` |
| FEM Não-Linear (validação analítica) | Yaw, L.L. — *3D Co-rotational Truss Formulation*, Walla Walla University (2011) | Fórmula fechada P(δ) para treliça rasa simétrica (Seção 10, Eq. 10.1), derivada do Crisfield — alvo de validação independente do FEM |

## A.2 `pso_rid.m` vs. Datta & Figueira (2011) — pontos confirmados e divergências

**Confirmado idêntico ao artigo:**
- Atualização de velocidade/posição real (Eq. 1–2).
- Contagem de dimensões (Eq. 5).
- Decodificação binário→inteiro (Eq. 6).
- Uso das regras de Deb para pbest/gbest (Seção 5, pág. 3628).
- Busca local no gbest (espírito da Seção 4.4).

**Divergências a resolver no Bloco 1:**

| # | Item | Artigo (Datta & Figueira) | Código atual | Ação proposta |
|---|---|---|---|---|
| D1 | Decodificação **discreta** | Índice direto 1..N via Eq. 6 (Seção 4.1) | Mapeamento proporcional customizado | Reverter para decodificação literal (índice direto), citando Seção 4.1/Eq. 6 no comentário; se o viés estatístico for comprovado pior, documentar com teste e decidir manter versão proporcional **como desvio citado**, não silencioso |
| D2 | Decodificação **inteira** | Bits pelo limite superior; fora de faixa tratado via violação de restrição, sem *clamping* no decoder | `min(decimal, config(k).max)`; ignora `var.min` | Remover o clamp interno; deixar o decodificador fiel à Eq. 6 e delegar limites ao tratamento de restrição (ou, no mínimo, corrigir para respeitar `min` também) |
| D3 | Velocidade **binária** (Tabelas 1–2, Eq. 8) | Assimétrica: `x=0 → v∈{0,1}`; `x=1 → v∈{0,-1}` | Simétrica `v∈{-1,0,1}` sempre + saturação pós-hoc | Reimplementar condicionando ao bit atual, citando Tabelas 1/2 linha a linha nos comentários |
| D4 | Mutação binária (`pm`) | Tie-break 50% dentro do subconjunto válido | Sorteio uniforme em `{-1,0,1}` (parte desperdiçada) | Ajustar para sortear apenas dentro do subconjunto válido do bit atual |
| D5 | Busca local em variável real | Mutação polinomial (Eq. 9), η auto-adaptativo em [25,45] | Gaussiana simples (5% do range) | Implementar Eq. 9 literalmente, citando o número da equação |

## A.3 Achados de validação do FEM Não-Linear (a documentar/testar no Bloco 2)

1. Peso da solução de referência Hadi recalculado: `2326.94 kg` vs `2325.2 kg` publicado
   (Δ≈0.075%) — não depende de carga/NR, só de geometria×área×densidade.
2. A solução de referência do Hadi, reavaliada pelo solver atual, viola o deslocamento
   (`51.02 mm > 50.80 mm`), de forma estável a variações de `n_inc` (10/50/200) — não é
   erro de convergência.
3. Achado extra desta sessão: caso de teste analítico fechado disponível (Yaw/Crisfield,
   Eq. 10.1) para uma treliça rasa simétrica — permite validar o Newton-Raphson contra uma
   solução exata, não apenas contra outro artigo.

## A.4 Itens fora de prioridade

- **Catálogo Awruch**: confirmado pelo usuário como catálogo de teste do professor, sem
  artigo de origem. Só limpeza pontual (remover `ref_areas`/`ref_peso` copiados do Hadi em
  `catalogo_awruch.m`), sem busca de referência.

---

# Parte B — Reorganização estrutural em 3 blocos

## B.1 Visão geral

```
CPIO III/
├── 01_pso_rid/                    # BLOCO 1
│   ├── pso_rid.m                  # solver PSO-RID fiel ao Datta & Figueira (2011)
│   └── (auxiliares só se necessário — preferir tudo dentro do próprio solver)
│
├── 02_fem_nao_linear/             # BLOCO 2 — dividido em 2 sub-blocos
│   ├── solver/
│   │   └── fem_nao_linear_solver.m   # recebe (nodes, elements, apoios, cargas, E, área...)
│   │                                  # devolve (weight, Sigma, u, diag) — SEM dados de nenhum problema específico embutido
│   └── problemas/
│       ├── problema_hadi_10barras.m      # só parâmetros (geometria/material/carga/catálogo) do caso Hadi
│       ├── problema_trelica_rasa_2barras.m # caso analítico Yaw/Crisfield (para teste unitário)
│       └── ... (novos problemas futuros = novos .m pequenos aqui)
│
├── 03_orquestrador/                # BLOCO 3
│   └── main_<nome_do_experimento>.m  # chama pso_rid + fem solver + problema, imprime relatório
│       (funções auxiliares só quando ocupam muitas linhas, ex.: blocos de texto/print)
│
├── 00_docs/
├── 04_resultados/
├── 05_legado/
│   └── pso_rid_pre_datta_2026-08-28.m   # cópia do pso_rid.m ANTES do ajuste de fidelidade
└── 06_testes/                      # testes matlab.unittest, um arquivo por bloco
```

Isso é muito próximo da estrutura atual (`01_src/pso`, `01_src/fem_truss`,
`02_casos_estudo`, `03_experimentos`) — a diferença principal é deixar explícito que:

1. O **solver FEM** (Bloco 2a) nunca conhece nenhum problema específico — só recebe
   parâmetros genéricos de treliça e devolve resultados. Isso já é quase verdade hoje
   (`fem_truss_nonlinear.m` tem a geometria do Hadi *hardcoded* dentro da função — **isso
   muda**: a geometria passa a ser parâmetro de entrada, não constante interna).
2. Os **problemas** (Bloco 2b) viram arquivos pequenos e independentes, um por caso de
   estudo, cada um só com dados (nós, elementos, apoios, cargas, catálogo/limites) — sem
   nenhuma lógica de solver.
3. O **orquestrador** (Bloco 3) concentra toda a lógica de "rodar o experimento" (loop de
   runs, chamar PSO, chamar FEM com o problema escolhido, montar relatório, plotar,
   salvar) em um único `.m` por experimento, evitando espalhar funções auxiliares em
   múltiplos arquivos — exceto blocos grandes de texto/impressão, que continuam isolados
   (ex.: `relatorio_comparativo.m`, `plot_convergencia.m`) para não poluir o fluxo
   principal de leitura.

## B.2 Citações teóricas embutidas no código

Cada fórmula relevante nos `.m` passa a levar um comentário citando **artigo + número da
equação/tabela**, não só o nome do artigo. Convenção proposta:

```matlab
% Eq. (14), Hadi & Alvani (2003): K_G = (P/L) * [mu^2 -lam*mu ...]
kg = (Pn / Ln) * [...];
```

```matlab
% Eq. (6), Datta & Figueira (2011): x = sum_{i=1}^{B} 2^(B-i) * b_i
decimal = round(bits) * pesos';
```

```matlab
% Regra de Deb (2000), pág. 316, critérios 1-3 (dominância feasível/violação)
if (pbest_viol(i) > 0) && (viol_atual < pbest_viol(i))
```

Isso vale para os 3 blocos, com prioridade para os pontos identificados na Parte A.2
(D1–D5) e para a montagem de K_E/K_G no solver FEM.

## B.3 Testes unitários por bloco (retomado de A, agora mapeado aos 3 blocos)

- **Bloco 1 (PSO-RID)**: `TestDecodificadorPSO.m`, `TestReglasDeb.m`,
  `TestVelocidadeBinaria.m` — testam `pso_rid.m` isolado, sem qualquer FEM.
- **Bloco 2 (FEM)**: `TestFemLinear.m`, `TestFemNaoLinear.m` — testam
  `fem_nao_linear_solver.m` isolado, alimentado por problemas sintéticos pequenos (não
  precisam ser um dos `problemas/*.m` reais) e pelo caso analítico fechado
  Yaw/Crisfield.
- **Bloco 3 (Orquestrador) / integração**: testes de ponta a ponta que chamam o
  orquestrador com um problema pequeno e conferem se o relatório/figura são gerados e se
  o melhor peso encontrado é plausível (ex.: dentro de X% da referência conhecida) —
  valida a integração dos 3 blocos juntos, não cada um isolado.

## B.4 Ordem de execução proposta

1. **Bloco 1 — PSO-RID** (mais autocontido, já com divergências D1–D5 mapeadas):
   - Copiar `pso_rid.m` atual para `05_legado/pso_rid_pre_datta_2026-08-28.m`.
   - Reescrever `pso_rid.m` corrigindo D1–D5, com citações inline (Seção B.2).
   - Validar rodando os experimentos existentes via Octave para conferir que ainda
     converge de forma razoável.
2. **Bloco 2 — FEM Não-Linear** (maior mudança estrutural: extrair geometria do Hadi de
   dentro da função para um arquivo de parâmetros):
   - Criar `fem_nao_linear_solver.m` genérico (recebe nós/elementos/apoios/cargas/E).
   - Criar `problema_hadi_10barras.m` com os dados atuais de `fem_truss_nonlinear.m`.
   - Criar `problema_trelica_rasa_2barras.m` para o teste analítico Yaw/Crisfield.
   - Repetir para a versão linear.
3. **Bloco 3 — Orquestrador**: reescrever os `run_*.m` como orquestradores únicos por
   experimento, consumindo os Blocos 1 e 2.
4. **Testes**: escrever os `matlab.unittest.TestCase` por bloco (Seção B.3), rodáveis
   isoladamente e em conjunto.

Este arquivo será atualizado conforme o plano evoluir; o `relatorio_mudancas_2026-08-28.md`
registra o que efetivamente foi executado, para comparação depois.
