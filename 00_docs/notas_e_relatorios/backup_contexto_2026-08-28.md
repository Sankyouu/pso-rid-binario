# Backup de Contexto — 2026-08-28

> Cópia de segurança da análise de validação teórica entregue na conversa de 2026-08-28,
> antes de iniciar a reorganização do projeto em blocos. Guardado "por desencargo de
> consciência" para não se perder caso o contexto da conversa seja resumido/perdido.

---

## 1. `pso_rid.m` vs. Datta & Figueira (2011) — checagem ponto a ponto

Artigo de referência: Dilip Datta, José Rui Figueira, *"A real-integer-discrete-coded
particle swarm optimization for design problems"*, Applied Soft Computing 11 (2011)
3625–3633. (`00_docs/artigos/Areal-integer-discrete-coded particle swarm optimization for
design problems.pdf`)

### Bate exatamente com o artigo

- Atualização de velocidade/posição real (Eq. 1–2) — idêntica.
- Contagem de dimensões (Eq. 5, `L = R + ΣB_j^I + ΣB_k^D`) — idêntica.
- Decodificação binário→inteiro (Eq. 6, `x = Σ 2^(B-i) b_i`) — o `pesos =
  2.^((n-1):-1:0)` do código reproduz essa soma exatamente.
- Uso das regras de Deb para atualizar pbest/gbest — o artigo (Seção 5, pág. 3628) diz
  textualmente que usa "the penalty-parameter-less constraint handling approach, proposed
  by Deb [8]" para tratar inviabilidade; o código faz o mesmo.
- Busca local no gbest (Seção 4.4: "a binary position of a particle is altered... with a
  small random probability") — presente em espírito no bloco "BUSCA LOCAL NO GBEST".

### Diverge do artigo (decisões silenciosas, não documentadas como desvio no código)

| Item | Datta & Figueira (2011) | `pso_rid.m` atual |
|---|---|---|
| Decodificação de variável **discreta** | Trata como inteiro-índice (1..N) via Eq. 6 direto (Seção 4.1) | Usa mapeamento **proporcional** customizado (`idx = floor(decimal*N/n_estados)+1`), descrito no comentário como correção de viés em relação à "versão anterior" — mudança deliberada, mas não é o que o artigo descreve |
| Decodificação de variável **inteira** | Bits dimensionados só pelo limite superior (igual ao código); fora de faixa parece ser tratado deixando o PSO empurrar via violação de restrição (ex.: `12 ≤ z ≤ 60` no problema da engrenagem, Eq. 12), não por *clamping* no decodificador | Código faz `min(decimal, config(k).max)` — clampa duro dentro do decodificador, e ignora `var.min` completamente (bug latente) |
| Atualização de velocidade **binária** (Tabelas 1 e 2, Eq. 8) | Regra **assimétrica**, condicionada ao bit atual: se `x=0`, `v ∈ {0,1}`; se `x=1`, `v ∈ {0,-1}` | Código usa regra **simétrica** `v ∈ {-1,0,1}` sempre, independente do bit atual, e depois satura `pos` em `[0,1]`. Parecido no resultado imediato, mas a memória de inércia (`vel(i,d)` usado no próximo `w*vel`) fica diferente da tabela do artigo em vários casos (ex.: linha 9 da Tabela 1, onde o artigo força `v=0` e o código pode guardar `v=-1`) |
| Mutação/exploração binária (`pm`) | Tie-break com 50% de probabilidade dentro do subconjunto válido ({0,1} ou {0,-1}) | Sorteio uniforme entre `{-1,0,1}` sem condicionar ao bit atual — parte dos sorteios é "desperdiçada" pelo clamping |
| Busca local em variável **real** | Mutação polinomial (Eq. 9) com índice de distribuição η auto-adaptativo em [25,45] | Perturbação **gaussiana** simples (`randn*sigma`, sigma=5% do range) |

Nenhum desses é necessariamente "errado" — são escolhas de engenharia razoáveis — mas como
Datta é a referência-norte, cada um precisa virar uma decisão explícita.

### Confirmação da regra de Deb na fonte primária

Fui direto na fonte primária — Kalyanmoy Deb, *"An efficient constraint handling method
for genetic algorithms"*, Comput. Methods Appl. Mech. Engrg. 186 (2000) 311–338 (já
presente em `00_docs/artigos/`), página 316:

> "1. Any feasible solution is preferred to any infeasible solution.
> 2. Among two feasible solutions, the one having better objective function value is
> preferred.
> 3. Among two infeasible solutions, the one having smaller constraint violation is
> preferred."

Isso bate **exatamente** com os comentários e a lógica condicional em `pso_rid.m` (regras
1/2/3 nas atualizações de pbest, gbest e aceitação da mutação). Validado contra a fonte
primária correta, não só contra o resumo do Datta.

---

## 2. Catálogo Awruch

Confirmado com o usuário: foi um catálogo montado pelo professor apenas para testar o
solver, sem artigo de origem — **não é prioridade de validação**. Fica só como limpeza
menor (tirar o `ref_areas`/`ref_peso` copiado do Hadi em `catalogo_awruch.m`, que é
enganoso porque nenhum desses valores pertence ao catálogo Awruch de fato).

---

## 3. Literatura adicional para embasar a FEM não-linear

Além do Crisfield Vol. 1 (já na pasta, citado tanto pelo Hadi quanto por outros), foi
encontrado:

**Louie L. Yaw, "3D Co-rotational Truss Formulation" (Walla Walla University, 2011)** —
nota técnica de livre acesso que deriva a mesma formulação corrotacional/Newton-Raphson do
Hadi (cita o próprio Crisfield como base), e traz um **exemplo com solução analítica
fechada**: uma treliça de 3 barras simétrica (tipo von Mises), com

```
P(δ) = (3EA/L³)·(z₀²δ + (3/2)z₀δ² + (1/2)δ³)
```

fórmula derivada no Crisfield (que já temos) para 1 barra e multiplicada por 3 pela
simetria. Isso é mais forte como validação teórica do que só bater com a Tabela 1 do
Hadi, porque dá um alvo fechado, calculável à mão, independente de qualquer solução de
outro algoritmo.

Fonte: <https://gab.wallawalla.edu/~louie.yaw/Co-rotational_docs/3Dcorot_truss.pdf>
(PDF baixado e lido localmente durante a sessão).

Proposta concreta (mantida no plano): montar um caso reduzido 2D análogo (treliça rasa de
2 barras, mesmo princípio geométrico) e usar essa fórmula fechada como oráculo em um teste
unitário do solver FEM não-linear.

---

## 4. Achados de validação já levantados (da rodada de leitura completa do projeto)

Cross-check de `fem_truss_nonlinear.m` / `catalogo_hadi.m` contra Hadi & Alvani (2003):

- Geometria, topologia, apoios, nós de carga: idênticos à Fig. 1 do artigo.
- E = 6.89×10⁴ MPa, densidade 2770 kg/m³, σ_max = 172.25 MPa, d_max = 50.80 mm: idênticos.
- Catálogo de 13 seções (65...29032 mm²): idêntico ao "segundo catálogo" do artigo.
- Matriz K_G do código bate termo a termo com a Eq. (14) do artigo.
- `caso.ref_areas` e `caso.ref_peso = 2325.2 kg` batem exatamente com a Tabela 1, Case 2
  discreto, do artigo.
- O refactor em `01_src` reproduz fielmente a lógica do código legado em `05_legado`
  (comparação linha a linha feita entre o antigo `truss.m` não-linear e o novo
  `fem_truss_nonlinear.m`).

Achados que precisam de investigação/documentação (evidência numérica, rodado via
Octave):

1. **Peso da solução de referência não bate exatamente**: `fem_truss_nonlinear(caso.ref_areas)`
   dá `weight = 2326.94 kg`, não `2325.2 kg` (Δ ≈ 0.075%). Como o peso só depende de
   área×comprimento×densidade (independente do carregamento/NR), a diferença é de
   geometria/arredondamento, não de convergência.
2. **A própria solução de referência do artigo viola o deslocamento nesta
   reimplementação**: `max|u| = 51.02 mm > d_max = 50.80 mm`. Testada sensibilidade a
   `n_inc` (10/50/200 — resultado idêntico) e à magnitude exata da carga (444822.16 N
   "exato" vs 445400 N "arredondado" do artigo — piora ainda mais). Não é bug de
   convergência do Newton-Raphson; é uma discrepância de formulação/arredondamento entre
   esta implementação e a do artigo, a isolar.
3. **`catalogo_awruch.m`** carrega `ref_areas`/`ref_peso` copiados do Hadi — confirmado
   numericamente que esses valores não pertencem às faixas do catálogo Awruch definido no
   mesmo arquivo. Código morto e enganoso (baixa prioridade, ver item 2 acima).
4. **Bug latente em `pso_rid.m`** no tipo `'I'` (inteiro): número de bits e decodificação
   usam só `var.max`, ignorando `var.min`. Hoje não afeta nada porque só o tipo `'D'` é
   usado nos casos de estudo, mas é uma bomba-relógio se `'I'` for usado no futuro.
5. **Melhor peso idêntico (2342.11 kg)** nas figuras salvas de linear e não-linear
   (`04_resultados/figuras/convergencia_hadi_*.png`) — pode ser coincidência (mesma
   semente, mesma combinação discreta próxima do ótimo) ou sinal de algo errado na
   separação dos casos. Precisa checar se os vetores de área são de fato iguais.
6. Não existe harness de teste automatizado — `debug_fem_nonlinear.m` é útil (rodei via
   Octave e confirma equilíbrio global, resíduo ~0.12 N) mas só imprime números, não
   falha/passa. `00_docs/notas_e_relatorios/` e `04_resultados/logs/` estavam vazios antes
   deste arquivo.

---

## 5. Estrutura de testes proposta (MATLAB `matlab.unittest.TestCase`)

Pasta de testes com uma classe por módulo, todas rodáveis via `runtests(...)`:

- **`TestFemLinear.m`** — equilíbrio `K*u == F`; caso estaticamente determinado (1 barra
  sob tração pura) com solução à mão; simetria geometria/carga → resultado espelhado.
- **`TestFemNaoLinear.m`** — equilíbrio com K_G; **caso da treliça rasa 3/2-barras vs.
  fórmula fechada de Yaw/Crisfield**; limite de carga pequena (NL → linear quando carga
  →0); invariância a `n_inc`; reprodução do benchmark Hadi Case 2 documentando o gap
  (peso 2326.94 vs 2325.2 kg; deslocamento 51.02 vs 50.80 mm) como teste de
  caracterização, não de igualdade estrita.
- **`TestDecodificadorPSO.m`** — `decodificar()` isolado contra Eq. 6 do Datta para
  padrões de bits conhecidos; viés estatístico da decodificação `'D'` por amostragem;
  exposição do bug do tipo `'I'`.
- **`TestReglasDeb.m`** — as 3 comparações do Deb (2000) isoladamente com pares
  (custo, violação) sintéticos.
- **`TestVelocidadeBinaria.m`** — reprodução de células das Tabelas 1 e 2 do Datta
  (combinações de `v_prev`, bit atual, pbest bit, gbest bit) vs. saída real do código.
- **`TestCustoTrelica.m`** — violação 0 dentro dos limites, crescimento quadrático fora,
  para valores sintéticos de `Sigma`/`u_livre`.
