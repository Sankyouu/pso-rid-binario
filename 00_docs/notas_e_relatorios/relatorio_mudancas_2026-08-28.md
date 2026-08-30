# Relatório de Mudanças — CPIO III

> Documento vivo, atualizado conforme cada mudança do `plano_final_reorganizacao_2026-08-28.md`
> é efetivamente aplicada. Cada linha registra estado anterior → mudança → justificativa
> (referência teórica) → status, para permitir comparação posterior (antes/depois).
>
> Status possíveis: `Planejado` · `Em andamento` · `Concluído` · `Revertido`

---

## Bloco 1 — PSO-RID

| # | Arquivo | Antes | Depois | Justificativa (referência) | Status |
|---|---|---|---|---|---|
| 1.1 | `01_src/pso/pso_rid.m` → cópia para `05_legado/pso_rid_pre_datta_2026_08_28.m` | — | Versão atual arquivada intacta (função renomeada para permitir chamada lado a lado) | Preservar estado pré-ajuste de fidelidade | **Concluído** |
| 1.2 | `pso_rid.m` — decodificação discreta (D1) | Mapeamento proporcional customizado `floor(x*N/2^B)+1` | Índice direto `x+1`, limites 1..N; overflow → inviável estrutural | [DF2011] Seção 4.1 / Eq. 6 | **Concluído** ⚠️ ver achado crítico |
| 1.3 | `pso_rid.m` — decodificação inteira (D2) | `min(decimal, max)`, ignora `min` | Valor direto da Eq. 6, sem clamp; `min` **e** `max` viram violação estrutural | [DF2011] Seção 4.1, Eq. 12 | **Concluído** |
| 1.4 | `pso_rid.m` — velocidade binária (D3) | Regra simétrica `v∈{-1,0,1}` + saturação | Regra assimétrica: `x=0 → v∈{0,+1}`; `x=1 → v∈{0,-1}` | [DF2011] Eq. 7–8, Tabelas 1–2 | **Concluído** |
| 1.5 | `pso_rid.m` — mutação binária (D4) | Sorteio uniforme `{-1,0,1}` | Com prob. `p_m` adota o valor alternativo dentro do subconjunto válido | [DF2011] Seção 4.3, col. 7 das Tab. 1–2 | **Concluído** |
| 1.6 | `pso_rid.m` — busca local (D5) | Gaussiana numa única dimensão sorteada | Percorre **todas** as dimensões: flip binário com prob. em ]0,0.05[ e mutação polinomial Eq. 9 nas reais | [DF2011] Seção 5 e Eq. 9 | **Concluído** |
| 1.7 | `pso_rid.m` — citações inline | Comentários genéricos | Cabeçalho com mapa equação→código; cada fórmula cita artigo + nº da eq./tabela/página | — | **Concluído** |
| 1.8 | `pso_rid.m` — dimensionamento de bits | `ceil(log2(limite+0.1))` para ambos os tipos | `I`: `ceil(log2(max+1))`; `D`: `ceil(log2(N))` — validado contra os exemplos numéricos do artigo (engrenagem 61→6 bits; mola 42→6, 70→7) | [DF2011] Seções 5.1 e 5.2 | **Concluído** |
| 1.9 | `pso_rid.m` — hiperparâmetros | `w` 0.9→0.4 linear, `c1=c2=2.0` | `w=1.0, c1=1.0, c2=2.0, pm=0.15` com auto-adaptação instantânea em `[0, valor]` (flag `auto_adaptativo`) | [DF2011] Seção 5, pág. 3628 | **Concluído** |
| 1.10 | `03_experimentos/run_hadi_nao_linear.m` | Passava `w_max`/`w_min` (não mais reconhecidos) | Atualizado para `w`/`c1`/`c2` com valores do artigo | Compatibilidade com 1.9 | **Concluído** |

**Resultado esperado / a preencher após execução:** comparação de peso ótimo, nº de
iterações até convergência e histórico de convergência entre a versão legado
(`05_legado/pso_rid_pre_datta_2026_08_28.m`) e a versão nova, rodando os mesmos
experimentos com a mesma semente.

### ⚠️ ACHADO CRÍTICO — consequência da decodificação literal (D1)

Ao implementar a Eq. 6 / Seção 4.1 do [DF2011] ao pé da letra (índice = valor binário + 1,
limites 1..N), surge uma consequência quantificada:

```
Catálogo Hadi: N = 13 opções  ->  B = ceil(log2(13)) = 4 bits  ->  2^4 = 16 estados
Códigos válidos por variável: 13/16 = 81.2%
Códigos inválidos (índice 14, 15, 16): 3/16 = 18.8%
Probabilidade de uma solução de 10 variáveis ser 100% válida: (13/16)^10 = 12.5%
```

Ou seja: **~87.5% das partículas geradas aleatoriamente são estruturalmente inviáveis**
(contêm ao menos um índice fora do catálogo) e não podem ser avaliadas pelo FEM.
Confirmado empiricamente no smoke test: 6 de 6 avaliações iniciais caíram em overflow.

**Este é exatamente o motivo pelo qual a versão anterior inventou o mapeamento
proporcional.** O mapeamento proporcional elimina 100% do overflow, ao custo de um viés
estatístico (opções do catálogo com probabilidades desiguais).

**Efeito colateral adicional a observar:** pela regra de [DEB2000] critério 3, entre duas
soluções inviáveis vence a de menor violação. A violação de overflow é um inteiro pequeno
(1 a ~30), enquanto a violação de tensão/deslocamento é enorme (~1e4 a 1e6 MPa²). Logo,
uma solução com índice apenas 1 unidade fora do catálogo **domina** uma solução
estruturalmente válida mas com tensão excessiva. Isso pode enviesar a busca. As duas
violações são de naturezas e escalas incomparáveis.

**Status da decisão:** implementada a versão literal (`'datta'`) como **padrão**, com o
mapeamento proporcional disponível como opção (`pso_params.decodificacao_discreta =
'proporcional'`) para comparação experimental. O overflow é medido e reportado em
`details.n_overflow`.

### Comparação empírica do D1 (isolada)

Experimento controlado: **mesmo código novo**, mesma semente (42), 100 partículas,
benchmark Hadi não-linear, **orçamento igualado em avaliações FEM** (~15 000) — a única
diferença entre A e B é a estratégia de decodificação:

| Config | Peso final | Viável? | Iterações | Avaliações FEM | Descartes por overflow | Tempo |
|---|---|---|---|---|---|---|
| A — `proporcional` | 2483.10 kg | sim | 150 | 15 150 | 0 | 450 s |
| B — `datta` (literal) | **2467.21 kg** | sim | 560 | 12 499 | 44 061 | 433 s |

**Conclusão preliminar: as duas estratégias são equivalentes na prática** (diferença de
0.6%, dentro do ruído de uma execução única). O desperdício de 78% das partículas no modo
`datta` é compensado pelo fato de que as avaliações descartadas são baratas (não chamam o
FEM), permitindo mais iterações no mesmo tempo de parede.

⚠️ **Importante — não é conclusão estatística.** Cada linha é **uma única execução**. Para
afirmar superioridade de qualquer das duas seria necessário um estudo multi-semente
(≥10 execuções por configuração) com média e desvio padrão. Isso fica registrado como
trabalho pendente.

### ⚠️ Achado secundário: possível regressão vinda de D3–D5/hiperparâmetros

Comparação adicional (semente 42, 100 partículas, 300 iterações, benchmark Hadi não-linear):

| Versão | Peso final | Iterações | Observação |
|---|---|---|---|
| Legado completo (pré-Datta) | **2356.80 kg** | 214 | `w` 0.9→0.4 linear, `c1=c2=2.0`, busca local gaussiana |
| Nova (fiel a Datta) | 2536.56 kg | 300 | `w,c1,c2,pm` auto-adaptativos em `[0, valor]` |

A versão legado obteve resultado **melhor**. Como o teste isolado do D1 mostrou que a
decodificação praticamente não influencia, a diferença deve vir de **D3–D5 e/ou dos
hiperparâmetros do artigo** (item 1.9). A hipótese mais provável é a auto-adaptação
instantânea: com `w, c1, c2 ~ U[0, valor_inicial]`, as médias efetivas caem para
`w≈0.5, c1≈0.5, c2≈1.0`, enfraquecendo bastante a atração ao pbest/gbest frente ao
`c1=c2=2.0` fixo da versão anterior.

**Como investigar (pendente):** rodar a versão nova com `auto_adaptativo = false` e
`w=0.9, c1=2.0, c2=2.0`, isolando o efeito dos hiperparâmetros do efeito das regras
D3–D5. Novamente, com múltiplas sementes.

**Nota metodológica:** ambas as comparações acima são de execução única e servem apenas
para levantar hipóteses, não para concluir.

---

## Bloco 2 — FEM Não-Linear (solver + parâmetros)

| # | Arquivo | Antes | Depois | Justificativa (referência) | Status |
|---|---|---|---|---|---|
| 2.1 | `02_fem_nao_linear/solver/fem_nao_linear_solver.m` | Geometria/material do Hadi *hardcoded* dentro de `fem_truss_nonlinear.m` | Solver genérico recebendo `problema` (nós, elementos, apoios, cargas, E, densidade) + `areas` + `opcoes` | Reuso para múltiplos problemas; [HA2003] Eq. 13–14–18 citadas no código | **Concluído** |
| 2.2 | `02_fem_nao_linear/solver/fem_linear_solver.m` | Idem em `fem_truss_linear.m` | Mesma separação solver/parâmetros | [HA2003] Sec. 4 | **Concluído** |
| 2.3 | `02_fem_nao_linear/problemas/problema_hadi_10barras.m` | — | Parâmetros do benchmark extraídos para arquivo próprio, com diagrama da Fig. 1 em ASCII e citação de cada valor | [HA2003] Fig. 1, Sec. 6.1, Tab. 1 | **Concluído** |
| 2.4 | `02_fem_nao_linear/problemas/problema_trelica_rasa_2barras.m` | — | Caso de validação analítica (treliça de von Mises) com solução fechada + carga crítica + inversor numérico | Derivação própria (ver abaixo) | **Concluído** |
| 2.5 | `01_src/fem_truss/fem_truss_*.m` | Solvers completos | **Cascas de compatibilidade** que chamam o solver genérico com o problema do Hadi | Manter `03_experimentos/` funcionando até o Bloco 3 | **Concluído** |
| 2.6 | `setup_paths.m` | Caminhos fixos | Inclui `02_fem_nao_linear/` e documenta os 3 blocos | — | **Concluído** |
| 2.7 | `02_casos_estudo/catalogo_awruch.m` | `ref_areas`/`ref_peso` copiados do Hadi (incorretos) | Removidos/marcados como "sem referência própria" | Confirmado: catálogo de teste do professor, sem artigo-fonte | Pendente (baixa prioridade) |

### Verificação de equivalência da migração

O refactor foi confirmado como **numericamente idêntico** ao código anterior — não é
aproximação, é diferença exatamente zero:

```
NAO LINEAR: dPeso=0.000e+00  dSigma=0.000e+00  du=0.000e+00
LINEAR    : dPeso=0.000e+00  dSigma=0.000e+00  du=0.000e+00
```

Há também um teste unitário permanente (`testCascaAntigaProduzResultadoIdentico`, em
ambas as classes de teste do FEM) que trava essa equivalência com `AbsTol = 0`.

### ✅ Validação analítica do Newton-Raphson — o resultado mais forte do Bloco 2

Foi construído um caso com **solução exata em forma fechada**: treliça rasa simétrica de
2 barras (treliça de von Mises). Para deslocamento vertical `w` do ápice:

```
L0   = sqrt(a² + h²)
L(w) = sqrt(a² + (h-w)²)
N(w) = E·A·(L(w) - L0)/L0            (deformação de engenharia)

P(w) = -2·N(w)·(h-w)/L(w)
```

Resultado da comparação (a=1000 mm, h=100 mm, E=70 GPa, A=100 mm²):

| Carga | `w` exato | `w` do solver | Erro relativo |
|---|---|---|---|
| 10% P_cr | 1.99269 mm | 1.99269 mm | 4.0e-12 |
| 25% P_cr | 5.23555 mm | 5.23555 mm | 8.4e-11 |
| 50% P_cr | 11.59273 mm | 11.59273 mm | 2.1e-12 |
| 75% P_cr | 20.14539 mm | 20.14539 mm | 7.5e-11 |
| 90% P_cr | 28.02480 mm | 28.02480 mm | 1.2e-11 |

Isso valida simultaneamente: atualização da geometria, matriz geométrica `K_G` (Eq. 14),
montagem das forças internas e convergência do Newton-Raphson — **contra um alvo
independente**, não contra o resultado de outro artigo.

Para provar que o caso é genuinamente não-linear (e não um teste trivial em que qualquer
solver acertaria), o mesmo teste confirma que a **análise linear erra 37.9%** nesse
regime (17.41 mm contra 28.02 mm exatos).

### ⚠️ Correção importante de método: qual fórmula analítica usar

O plano original previa usar a fórmula de Yaw (2011) Eq. 10.1 / Crisfield como oráculo:

```
P(δ) = (3EA/L³)·(z₀²δ + (3/2)z₀δ² + (1/2)δ³)
```

**Isso teria sido um erro metodológico.** Essa fórmula clássica é derivada com
**deformação de Green** e hipótese de abatimento, enquanto o solver deste projeto usa
**deformação de engenharia** `(L-L0)/L0`. As duas não coincidem fora do regime de
deformações infinitesimais — o desacordo resultante seria erroneamente interpretado como
bug no solver.

Por isso a fórmula acima foi **derivada do zero para a mesma medida de deformação do
código**, a partir do equilíbrio vertical do nó do ápice na configuração deformada. Essa
observação está registrada como comentário tanto no solver quanto no arquivo do problema,
para não se perder.

### Nota sobre o ponto limite (snap-through)

A treliça rasa tem ponto limite: `P(w)` cresce, atinge um máximo (`P_crítica`) e decresce
até zero em `w = h`. Um solver de **controle de carga** como o do Bloco 2a só consegue
seguir o ramo estável abaixo dessa carga. O arquivo do problema calcula `P_critica` e o
inversor `w_de_P` **recusa explicitamente** cargas acima dela (com erro nomeado), em vez
de devolver um número silenciosamente errado. Há teste para isso.

### Achados anteriores — reconfirmados após a migração

Ambos permanecem exatamente iguais (como esperado, por ser refactor estrutural):
- Peso da solução de referência Hadi: **2326.94 kg** vs **2325.2 kg** publicado (+0.075%).
- Deslocamento máximo: **51.02 mm** vs limite **50.80 mm** (+0.43%).

Estão travados como **testes de caracterização** (`testHadi_CaracterizacaoDaSolucaoDeReferencia`),
que falham se a discrepância mudar de tamanho — documentando o desvio em vez de escondê-lo.

---

## Bloco 3 — Orquestrador

| # | Arquivo | Antes | Depois | Justificativa | Status |
|---|---|---|---|---|---|
| 3.1 | `03_orquestrador/main_hadi_nao_linear.m` | `run_hadi_nao_linear.m` chamava `custo_trelica.m` externo | Orquestrador autocontido: função objetivo, laço multi-start, seleção do melhor e estatísticas como funções locais | Diretriz: um `.m` por experimento, tudo embutido exceto blocos longos de saída | **Concluído** |
| 3.2 | `03_orquestrador/main_hadi_linear.m` | idem | idem, + confronto do projeto ótimo linear reavaliado em análise não linear | Materializa o argumento de [HA2003] Sec. 1 | **Concluído** |
| 3.3 | `03_orquestrador/main_awruch_discreto.m` | idem | idem, + verificação de que cada área pertence ao catálogo da sua própria barra | Catálogos de cardinalidades diferentes exercitam a Eq. (5) por variável | **Concluído** |
| 3.4 | `03_orquestrador/main_estudo_estatistico.m` | `run_estudo_particulas.m` (50 vs 100 partículas, execuções únicas) | **Novo**: comparação multi-semente pareada de 4 configurações, com orçamento igualado em avaliações FEM | Responder Q1 e Q2 deixadas em aberto no Bloco 1 | **Concluído** |
| 3.5 | `03_orquestrador/auxiliares/` | `01_src/utils/` | `relatorio_comparativo.m`, `plot_convergencia.m`, `salvar_figura.m` movidos para cá | São exatamente os "blocos longos de impressão/gráfico" que a diretriz manda manter fora | **Concluído** |
| 3.6 | `01_src/utils/custo_trelica.m` | Função externa compartilhada | **Removida**; embutida como `avaliar_projeto` local em cada orquestrador | Diretriz de orquestrador autocontido | **Concluído** |

### Reestruturação final de pastas

A renumeração adiada no Bloco 2 foi concluída (`git mv`, histórico preservado):

| Antes | Depois |
|---|---|
| `01_src/pso/` | `01_pso_rid/` |
| `01_src/fem_truss/` | removido (cascas de compatibilidade não são mais necessárias) |
| `01_src/utils/` | `03_orquestrador/auxiliares/` |
| `02_casos_estudo/` | absorvido em `02_fem_nao_linear/problemas/` |
| `03_experimentos/` | `03_orquestrador/` |

Arquivados em `05_legado/` para comparação: `pso_rid_pre_datta_2026_08_28.m`,
`fem_truss_nonlinear_pre_bloco2.m`, `fem_truss_linear_pre_bloco2.m` e
`pre_bloco3/` (os `run_*.m`, `debug_fem_nonlinear.m` e os catálogos antigos).

Os solvers FEM originais foram recuperados do histórico do git e sua equivalência com os
solvers genéricos permanece travada por teste (`AbsTol = 0`) — agora comparando contra o
**código original arquivado**, não contra uma casca.

### Decisão sobre a duplicação da função objetivo

`avaliar_projeto` aparece em cada orquestrador (~15 linhas). Isso é duplicação
deliberada, seguindo a diretriz de orquestrador autocontido: o custo de ler quatro cópias
pequenas é menor que o de perseguir uma função compartilhada externa para entender um
experimento. A fórmula é idêntica nos quatro e cita [HA2003] Eq. (1) e (2).

### O estudo estatístico (item que estava pendente do Bloco 1)

`main_estudo_estatistico.m` compara 4 configurações isolando **uma variável por vez**:

| Config | Decodificação | Hiperparâmetros | Responde |
|---|---|---|---|
| A | `datta` | auto-adaptativos [DF2011] | (base) |
| B | `proporcional` | auto-adaptativos | **Q1** (efeito da decodificação) |
| C | `datta` | fixos (`w=0.9, c1=c2=2.0`) | **Q2** (efeito da auto-adaptação) |
| D | `proporcional` | fixos | comportamento legado |

Metodologia: sementes **pareadas** (todas as configurações veem as mesmas), orçamento
igualado em **avaliações FEM** (não em iterações — configurações com muito descarte
estrutural fariam mais iterações no mesmo tempo, enviesando a comparação), e leitura
descritiva comparando a diferença observada contra o desvio padrão. O relatório indica
explicitamente quando a diferença é menor que o desvio, em vez de sugerir conclusão
onde não há evidência. Logs salvos em `04_resultados/logs/`.

#### Correção do mecanismo de parada por orçamento

A primeira implementação devolvia **custo infinito** ao esgotar o orçamento, deixando a
parada por conta do critério de estagnação. Isso era implícito e frágil: bastava afrouxar
`tol_estagnacao` para o laço girar dezenas de milhares de iterações sem realizar nenhuma
avaliação útil — foi exatamente o que aconteceu num teste de fumaça com
`tol_estagnacao = 1e9`.

Substituído por parada **explícita**: ao esgotar o orçamento a função objetivo lança
`main_estudo_estatistico:orcamentoEsgotado`, capturado pelo orquestrador. Como o erro
interrompe o `pso_rid` antes do retorno normal, a closure mantém por conta própria o
melhor resultado já visto, aplicando a mesma regra de dominância de Deb do solver
(`rid_domina_deb`) — assim nada se perde ao abortar.

Coberto por `TestIntegracao/testEstudoEstatistico_RespeitaOrcamentoDeAvaliacoes`, que
usa `tol_estagnacao = 1e9` de propósito: se a parada voltar a depender da estagnação, o
teste falha por tempo.

#### Correção metodológica: análise pareada, não comparação de desvios marginais

A primeira versão de `imprimir_conclusoes` comparava a diferença das médias contra os
**desvios marginais** de cada configuração. **Esse é o teste errado para amostras
pareadas.** Os desvios marginais carregam a variabilidade *entre sementes* — exatamente a
componente que o pareamento existe para eliminar.

O correto é trabalhar com o vetor de diferenças por semente
`d_i = peso_cfg2(i) − peso_cfg1(i)` e comparar `média(d)` com o **erro padrão da média**,
`desvio(d)/√n` (estatística t pareada), complementado pela contagem de vitórias, que é
robusta a valores extremos. Implementado com tabela de valores críticos embutida, para
não depender do Statistics Toolbox.

---

# ✅ RESULTADO DO ESTUDO ESTATÍSTICO — Q1 e Q2 respondidas

Execução: **4 configurações × 10 sementes pareadas**, orçamento de 15.000 avaliações FEM,
análise **linear**, em 129 s. Log: `04_resultados/logs/estudo_estatistico_linear_20260828_145632.txt`

| Config | Melhor | Média | Desvio | Viabilidade |
|---|---|---|---|---|
| A — Datta puro | 2374,81 | 2551,24 | 173,54 | 100% |
| B — A + decodificação proporcional | **2342,11** | 2467,37 | **88,78** | 100% |
| C — A + hiperparâmetros fixos | 2390,64 | 2582,54 | 230,14 | 100% |
| D — fixos + proporcional (legado) | 2380,44 | 2533,60 | 120,84 | 100% |

Referência `[HA2003]`: 2325,2 kg.

### Análise pareada

| Comparação | Média das diferenças | Erro padrão | t pareado | Vitórias | Conclusão |
|---|---|---|---|---|---|
| **Q1** B − A (decodificação) | −83,87 kg | 65,82 | −1,274 | 7/10 | **Sem evidência** (\|t\| < 2,262) |
| **Q2** C − A (hiperparâmetros) | +31,30 kg | 49,57 | +0,631 | 4/10 | **Sem evidência** |
| Extra D − A (legado) | −17,64 kg | 57,21 | −0,308 | 5/10 | **Sem evidência** |

### Conclusões

**Q1 — a decodificação discreta não importa estatisticamente neste problema.** Apesar de
o modo `datta` descartar ~87% das partículas, isso é compensado por os descartes serem
baratos (não chamam o FEM). A tendência favorece o proporcional (7 vitórias em 10), mas
não atinge significância.

**Q2 — a hipótese levantada no Bloco 1 está REFUTADA.** Eu havia suspeitado que a
auto-adaptação instantânea dos hiperparâmetros (`[DF2011]` Sec. 5) prejudicava a
convergência, com base na observação de execução única (legado 2356,80 kg vs nova
2536,56 kg). Com 10 sementes pareadas, a diferença é de +31 kg com t = 0,63 — **puro
ruído**. A observação original era um artefato de amostra de tamanho 1.

**Observação prática que sobrevive:** a decodificação proporcional é sensivelmente **mais
consistente** (desvio 88,78 contra 173,54 do `datta`). Média estatisticamente
indistinguível, mas menor dispersão — o que tem valor prático mesmo sem significância na
média.

**Recomendação:** manter `'datta'` como padrão (fidelidade ao artigo-norte, sem custo
demonstrável) e registrar a menor dispersão do `'proporcional'` como alternativa
justificada quando reprodutibilidade importar mais que fidelidade literal.

### Ressalvas honestas

- Estudo feito com FEM **linear**. As conclusões podem não transferir para o não-linear
  (mais caro: ~4 h para o mesmo desenho). Pendente.
- n = 10 detecta apenas efeitos grandes. Um efeito real de ~50 kg passaria despercebido.
- O melhor resultado global (2342,11 kg, config B) coincide **exatamente** com o valor
  registrado nas figuras salvas do código legado — verificação de consistência
  independente entre a implementação nova e a antiga.

---

## Testes de integração (Bloco 3)

`06_testes/TestIntegracao.m` cobre o que só aparece na junção dos blocos:

- interface do arquivo de problema servindo simultaneamente Bloco 1 e Bloco 2;
- toda área decodificada pelo PSO é aceita pelo FEM (ciclo fechado);
- solução final pertence ao catálogo (e, no caso Awruch, ao catálogo **da barra certa**);
- dimensionamento independente de bits para catálogos de cardinalidades diferentes;
- regressão: `problema_awruch_10barras` **não** herdou a referência incorreta do Hadi;
- mesmo struct de problema aceito pelos dois solvers sem adaptação;
- peso devolvido pelo PSO confere com o recalculado pelo FEM;
- histórico do g-best não crescente;
- **reprodutibilidade por semente** (requisito para o estudo pareado fazer sentido).

---

## Testes unitários (novo — não existia antes)

Framework: `matlab.unittest` nativo, conforme solicitado. **Nota de ambiente:** o Octave
8.4 (também instalado) **não** implementa `matlab.unittest` — os testes exigem MATLAB.
O Octave continua útil para rodar os solvers diretamente.

Executar tudo: `matlab -batch "cd('06_testes'); run_todos_testes"` (ou `runtests('06_testes')`)

**Status atual da suite: 53 testes, 53 passando, 0 falhas — 4.73 s.**

| # | Arquivo | Bloco | O que cobre | Status |
|---|---|---|---|---|
| T1 | `06_testes/TestDecodificadorPSO.m` | 1 | Eq. 5 (dimensionamento de bits, conferido contra os exemplos numéricos das Seções 5.1/5.2 do artigo), Eq. 6 (soma ponderada), D1 (ambos os modos + caracterização do trade-off), D2 (limites inferior **e** superior) | **Concluído — 16/16 passam** |
| T2 | `06_testes/TestRegrasDeb.m` | 1 | Os 3 critérios de [DEB2000] pág. 316 isoladamente, inicialização com Inf/Inf, e propriedades estruturais (antissimetria, irreflexividade, transitividade exaustiva) | **Concluído** |
| T3 | `06_testes/TestVelocidadeBinaria.m` | 1 | Eq. 7 (domínio admissível por bit — regressão do D3), Eq. 8 (mapeamento determinístico), coluna 7 das Tabelas 1–2 (frequência de `p_m` verificada estatisticamente), empate 50/50 | **Concluído** |
| T4 | `06_testes/TestMutacaoPolinomial.m` | 1 | Eq. 9: limites do domínio, propriedade `η` maior → perturbação menor, localidade na faixa [25,45] do artigo, centragem, casos degenerados | **Concluído** |
| T5 | `06_testes/TestFemLinear.m` | 2 | — | Planejado |
| T6 | `06_testes/TestFemNaoLinear.m` | 2 (inclui caso analítico Yaw/Crisfield) | — | Planejado |
| T7 | `06_testes/TestCustoTrelica.m` | 2/3 | — | Planejado |
| T8 | `06_testes/TestIntegracao.m` | 1+2+3 | — | Planejado |

### Refatoração de testabilidade (não estava no plano original)

Para atender ao requisito de *"analisar cada pedaço do solver isoladamente"*, os
auxiliares do PSO foram **extraídos de funções locais para arquivos próprios**, ficando
diretamente chamáveis pelos testes (antes só seriam alcançáveis por indireção):

| Arquivo novo | Equação/referência | Substitui |
|---|---|---|
| `01_src/pso/rid_mapear_dimensoes.m` | Eq. (5) [DF2011] | bloco inline em `pso_rid.m` |
| `01_src/pso/rid_decodificar.m` | Eq. (6) + Sec. 4.1 [DF2011] | função local `decodificar` |
| `01_src/pso/rid_velocidade_binaria.m` | Eq. (7)–(8), Tab. 1–2 [DF2011] | bloco inline em `pso_rid.m` |
| `01_src/pso/rid_mutacao_polinomial.m` | Eq. (9) [DF2011] | função local `mutacao_polinomial` |
| `01_src/pso/rid_domina_deb.m` | pág. 316 [DEB2000] | função local `domina_deb` |

`pso_rid.m` ficou reduzido ao **loop do PSO**, chamando esses auxiliares. Isso é coerente
com a decisão do Bloco 3 (orquestrador concentrado), mas aplicada ao Bloco 1: o solver
concentra o fluxo, e cada fórmula citável vive num arquivo pequeno e testável.

---

## Infraestrutura

| # | Mudança | Justificativa | Status |
|---|---|---|---|
| I1 | `git init` + commit baseline (`d0c0302`) | Projeto não tinha controle de versão; necessário para reorganização segura e reversível | **Concluído** (2026-08-28) |
