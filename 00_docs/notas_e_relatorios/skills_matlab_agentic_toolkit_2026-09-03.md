# Skills do MATLAB Agentic Toolkit aplicáveis ao projeto CPIO III

Levantamento feito em 2026-09-03, após leitura completa do projeto (Blocos 1-3 e
`06_testes/`). Lista as skills do toolkit que podem melhorar o projeto **respeitando
as diretrizes já estabelecidas** (sem over-engineering, sem alterar `05_legado/`,
preservando os identificadores de erro usados pelos testes, sem substituir o
PSO-RID por solvers prontos).

---

## Alta relevância — mapeiam para itens já em aberto no código

| Skill | Por quê |
|---|---|
| `matlab-software-development:matlab-optimize-performance` | `01_pso_rid/pso_rid.m` ainda tem laços aninhados não vetorizados (`for i=1:n_particulas` / `for d=1:total_dim`, linhas ~290-330 e ~435-497), enquanto `02_fem_nao_linear/fem_nao_linear_solver.m` já passou por esse processo (comentário "profile de 2026-09-02"). O workflow baseline→profile→otimizar→verificar é o mesmo que o projeto já aplicou ao FEM. |
| `parallel-computing:matlab-diagnose-parfor` | `03_orquestrador/main_estudo_estatistico.m` tem lógica própria de limite de workers por memória (`resolver_workers`, `preparar_pool`) — candidato natural a diagnóstico/simplificação com a skill dedicada. |
| `matlab-data-import-and-analysis:matlab-analyze-data` | O próprio código deixa isso como pendência: *"Para um teste livre de distribuição, aplique Wilcoxon pareado (signrank)"* (`main_estudo_estatistico.m`). Statistics and Machine Learning Toolbox está instalado — fecha um TODO já declarado, não é scope creep. |
| `matlab-programming:matlab-validate-function-arguments` | `pso_rid.m`, `fem_linear_solver.m`, `fem_nao_linear_solver.m` validam entradas manualmente com `assert`/identificadores de erro customizados que os testes checam (`verifyError(..., 'fem_nao_linear_solver:areasIncompativeis')`). Migrar para blocos `arguments` exige preservar esses identificadores — a skill cobre isso. |
| `matlab-software-development:matlab-write-performance-tests` | Não há guarda de regressão para o ganho de desempenho já documentado na vetorização do FEM (`accumarray`, máscara lógica em vez de `setdiff`). Complementa `06_testes/` no mesmo espírito de "trava numérica" que o projeto já usa. |
| `matlab-core:matlab-debug-code` | Há divergências abertas e não resolvidas: taxa de sucesso de 1% nas engrenagens vs. 100% do artigo, D6 (reinicialização) nunca medido no FEM não linear. Workflow estruturado de debug ajuda a isolar causa sem repetir hipóteses já descartadas. |

## Média relevância — boas práticas, sem violar as diretrizes do projeto

| Skill | Por quê |
|---|---|
| `matlab-core:matlab-review-code` | Passagem de revisão sistemática — compatível com a filosofia de "sem abstração prematura" que o próprio projeto já pratica. |
| `matlab-software-development:matlab-modernize-code` | Aplicável ao código ativo (pso_rid, solvers, orquestradores) — não a `05_legado/`, que é mantido congelado de propósito para comparação. |
| `math-and-optimization:matlab-use-symbolic-math` | Poderia verificar simbolicamente derivações fechadas já usadas como oráculo (fórmula da treliça rasa de von Mises, ou a redundância g7≡0 da mola) — reforça a cultura de "caracterização com evidência" que os testes já seguem. |
| `matlab-core:matlab-write-help` | Formalizar H1-line/help para que `help pso_rid` funcione — os comentários já são extensos, só faltaria o formato padrão. |

## Baixa relevância — só se o usuário quiser expandir escopo

| Skill | Por quê |
|---|---|
| `reporting-and-database-access:matlab-generate-report` | Substituiria os blocos `fprintf` manuais de `relatorio_comparativo.m` por PDF/HTML — mudança estrutural maior, não pedida. |
| `matlab-software-development:matlab-package-toolbox` | Só faz sentido se houver intenção de distribuir o projeto como toolbox instalável. |
| `matlab-core:matlab-create-live-script` | Documentação/demo interativa — o formato "plain text live code" exige R2025a+, e a instalação aqui é R2024a. |

## Explicitamente NÃO recomendadas — conflitam com o propósito do projeto

- **`math-and-optimization:matlab-solve-optimization` / Global Optimization Toolbox** — usar solvers prontos do MATLAB contrariaria o objetivo central: o projeto existe para implementar e validar fielmente o PSO-RID de `[DF2011]` contra a literatura, não para substituí-lo.
- **`model-based-design-core:*`, `model-based-system-engineering:*`, `simulink-*`, `verification-validation-and-test:*`, `matlab-app-building:*`** — não há nenhum modelo Simulink nem GUI no projeto.
- **`code-generation:*`, `ai-and-statistics:matlab-train-network/classify-tabular-data/etc.`, e todos os toolkits de domínio** (aerospace, automotive, RF, wireless, radar, imagem, biologia/finanças computacional) — fora do escopo de otimização estrutural de treliças.
