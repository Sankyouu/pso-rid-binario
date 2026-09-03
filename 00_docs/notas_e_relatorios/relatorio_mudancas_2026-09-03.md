# Vetorização, paralelismo opcional e validação de argumentos — 2026-09-03

Aplicado a `pso_rid.m`, aos dois solvers FEM e aos 7 orquestradores com caso de
estudo. `main_estudo_estatistico.m`, `06_testes/` e `05_legado/` ficaram fora
do escopo. Motivação completa de cada mudança (por que cada skill do
matlab-agentic-toolkit se aplica) em
`skills_matlab_agentic_toolkit_2026-09-03.md`, nesta mesma pasta.

## O que mudou

**Solvers FEM**
- `fem_linear_solver.m`: montagem da matriz de rigidez vetorizada, reusando o
  padrão `accumarray`/índices lineares já usado em `fem_nao_linear_solver.m`.
- Os dois solvers ganharam bloco `arguments` para checar tipo e forma da
  entrada.

**pso_rid.m**
- Seção C (atualização de velocidade/posição) vetorizada por completo.
- Seção A (avaliação de partículas) separada em `parfor` opcional
  (`pso_params.n_workers`, padrão 0 = serial). Não usa números aleatórios,
  então ligar paralelismo aqui não afeta reprodutibilidade.
- Bloco `arguments`, preservando a chamada-sentinela `pso_rid('auxiliares')`.

**Orquestradores**
- Laço multi-start (`for r = 1:n_runs`) paralelizável via novo argumento
  posicional `n_workers` (padrão 0 = serial).
- `garantir_caminhos()`, antes duplicada em cada arquivo, e a gestão de pool
  paralelo (extraída de `main_estudo_estatistico.m`) viraram helpers
  compartilhados: `03_orquestrador/auxiliares/garantir_caminhos.m` e
  `preparar_pool.m`.
- Bloco `arguments` substitui os `if nargin < N ...` manuais.

## Consequências aceitas

- **Reprodutibilidade por semente mudou.** A vetorização da Seção C do
  `pso_rid` e a semente própria por run (necessária para o laço multi-start
  ser seguro em paralelo) alteram a sequência de números aleatórios
  consumida. A mesma seed não reproduz mais bit a bit resultados gerados
  antes desta data — inclusive os já publicados em `04_resultados/`.
- Identificadores de erro de validação de `modo` nos orquestradores (ex.:
  `main_hadi_20barras:modoInvalido`) viraram genéricos do MATLAB
  (`mustBeMember`); nenhum teste dependia deles.
- `preparar_pool.m` limita workers pela memória disponível, não pelo número
  de núcleos — mesmo critério já usado em `main_estudo_estatistico.m`, para
  não travar a máquina em execução paralela.

## Evidência

Suíte completa: 109/110 testes passam. A única falha é
`TestFemLinear/testSolverOriginalProduzResultadoIdentico` — comparação
`AbsTol 0` contra o solver legado, quebrada pela reordenação de soma via
`accumarray` (divergência ~1e-14). É o mesmo precedente já documentado em
`TestFemNaoLinear.m` quando aquele solver foi vetorizado.

Rodados end-to-end sem erro: `main_hadi_linear`, `main_hadi_nao_linear` e
`main_datta_engrenagens` (3 runs paralelas por semente, ótimo global
reencontrado).
