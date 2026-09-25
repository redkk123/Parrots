# Longevidade em Psittaciformes: história evolutiva e bases genômicas

Protocolo completo: [`docs/protocol.md`](docs/protocol.md) (versão de 25 set. 2026).

**Estado atual:** Fase 0 (auditoria de novidade e viabilidade) executada parcialmente. Nenhuma análise molecular foi feita. **Decisão provisória: cenário B** (ver [`results/phase0/phase0_feasibility.md`](results/phase0/phase0_feasibility.md)).

## Entregáveis do primeiro ciclo

| Entregável do protocolo | Arquivo | Estado |
| --- | --- | --- |
| `literature.tsv` | [`data/curated/literature.tsv`](data/curated/literature.tsv) | 27 referências; maioria verificada só por resumo |
| `novelty_matrix.md` | [`results/phase0/novelty_matrix.md`](results/phase0/novelty_matrix.md) | preliminar (sem busca sistemática) |
| `species.tsv` | [`results/phase0/species.tsv`](results/phase0/species.tsv) | 413 táxons da *supertree* |
| `phenotypes.tsv` | [`results/phase0/phenotypes.tsv`](results/phase0/phenotypes.tsv) | 218 spp. com expectativa de vida; longevidade máxima pendente |
| `assemblies.tsv` | [`results/phase0/assemblies.tsv`](results/phase0/assemblies.tsv) | VGP verificado; NCBI pendente |
| `cancer_comparisons.tsv` | [`data/curated/cancer_comparisons.tsv`](data/curated/cancer_comparisons.tsv) | contrastes publicados + plano por filogenia |
| `neurodegeneration_evidence_map.tsv` | [`data/curated/neurodegeneration_evidence_map.tsv`](data/curated/neurodegeneration_evidence_map.tsv) | matriz de lacunas |
| `oxidative_resilience.tsv` | [`data/curated/oxidative_resilience.tsv`](data/curated/oxidative_resilience.tsv) | por etapa (produção → tolerância) |
| `integrated_axes_summary.md` | [`results/phase0/integrated_axes_summary.md`](results/phase0/integrated_axes_summary.md) | — |
| `species_tree.nwk` | [`results/phase0/species_tree.nwk`](results/phase0/species_tree.nwk) (218 spp.) e [`supertree_burgio2019_full.nwk`](results/phase0/supertree_burgio2019_full.nwk) | árvore alternativa (Smith et al.) pendente |
| `phase0_feasibility.md` | [`results/phase0/phase0_feasibility.md`](results/phase0/phase0_feasibility.md) | — |
| `analysis_plan.md` | [`results/phase0/analysis_plan.md`](results/phase0/analysis_plan.md) | v1.0 |

Resultados auxiliares em `results/phase0/`: ajustes alométricos (`model_fits.tsv`), sensibilidades, estados ancestrais, vizinhos filogenéticos de *Melopsittacus*, simulação de poder e figuras.

## Reproduzir

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/00_fetch_external.sh              # clona fontes em commits fixos (config/sources.tsv)
.venv/bin/python scripts/01_build_phase0_tables.py
.venv/bin/python scripts/02_phase0_analysis.py
```

R não estava disponível neste ambiente, então as análises filogenéticas (PGLS com λ de Pagel/OU, reconstrução ancestral BM) foram implementadas em Python (`scripts/lib/phylo.py`). Antes do artigo, recomenda-se validá-las contra `ape`/`phytools`/`nlme`.

## Estrutura

```
config/     fontes externas e versões fixadas
data/
  curated/  tabelas curadas à mão (com nível de verificação por linha)
  external/ downloads (não versionado)
docs/       protocolo
scripts/    pipeline numerado + lib/
results/    saídas por fase
logs/       logs e checksums das fontes
manuscript/ (vazio)
```

## Dados de terceiros

Os valores por espécie de expectativa de vida, massa e cérebro são **derivados** dos resultados processados de Smeele et al. (2022; [repositório](https://github.com/simeonqs/Coevolution_of_relative_brain_size_and_life_expectancy_in_parrots)). Os registros ZIMS brutos não são públicos. A árvore é a de Burgio et al. (2019), obtida via [`parrottraits`](https://github.com/trashbirdecology/parrottraits). Cite as fontes originais.

## Limitações deste ciclo

- A política de rede bloqueou NCBI, Crossref, Europe PMC, Dryad, Zenodo, Figshare, bioRxiv e CRAN. Por isso ficaram pendentes: accessions e BUSCO, o mapeamento dos 94 + 22 genomas short-read, a árvore de Smith et al. e a busca bibliográfica sistemática.
- Linhas marcadas com `to_verify` ou `unverified` **não** devem ser usadas como fato.
