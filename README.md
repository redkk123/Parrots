# Longevidade em Psittaciformes: história evolutiva e bases genômicas

Protocolo completo: [`docs/protocol.md`](docs/protocol.md) (versão de 25 set. 2026).

**Estado atual:** Fase 0 e Fase 0b (mapeamento genoma × fenótipo, desfecho alternativo, mudanças de regime) concluídas. Nenhuma análise molecular foi feita ainda. **Decisão: cenário B, agora com amostra viável (105 espécies)**. Ver [`results/phase0/phase0_feasibility.md`](results/phase0/phase0_feasibility.md) e [`results/phase0b/phase0b_report.md`](results/phase0b/phase0b_report.md).

## Entregáveis do primeiro ciclo

| Entregável do protocolo | Arquivo | Estado |
| --- | --- | --- |
| `literature.tsv` | [`data/curated/literature.tsv`](data/curated/literature.tsv) | 31 referências/fontes; maioria verificada só por resumo |
| `novelty_matrix.md` | [`results/phase0/novelty_matrix.md`](results/phase0/novelty_matrix.md) | preliminar (sem busca sistemática) |
| `species.tsv` | [`results/phase0/species.tsv`](results/phase0/species.tsv) | 413 táxons da *supertree* |
| `phenotypes.tsv` | [`results/phase0/phenotypes.tsv`](results/phase0/phenotypes.tsv) + [`results/phase0b/species_data_matrix_with_residuals.tsv`](results/phase0b/species_data_matrix_with_residuals.tsv) | 218 spp. com expectativa de vida; 173 com longevidade máxima (AnAge/Amniote) |
| `assemblies.tsv` | [`results/phase0/assemblies.tsv`](results/phase0/assemblies.tsv) + [`results/phase0b/species_data_matrix.tsv`](results/phase0b/species_data_matrix.tsv) | VGP verificado; dados WGS por espécie via metadados SRA; accessions de montagem NCBI pendentes |
| `cancer_comparisons.tsv` | [`data/curated/cancer_comparisons.tsv`](data/curated/cancer_comparisons.tsv) | contrastes publicados + plano por filogenia |
| `neurodegeneration_evidence_map.tsv` | [`data/curated/neurodegeneration_evidence_map.tsv`](data/curated/neurodegeneration_evidence_map.tsv) | matriz de lacunas |
| `oxidative_resilience.tsv` | [`data/curated/oxidative_resilience.tsv`](data/curated/oxidative_resilience.tsv) | por etapa (produção → tolerância) |
| `integrated_axes_summary.md` | [`results/phase0/integrated_axes_summary.md`](results/phase0/integrated_axes_summary.md) | — |
| `species_tree.nwk` | [`results/phase0/species_tree.nwk`](results/phase0/species_tree.nwk) (218 spp.) e [`supertree_burgio2019_full.nwk`](results/phase0/supertree_burgio2019_full.nwk) | árvore alternativa (Smith et al.) pendente |
| `phase0_feasibility.md` | [`results/phase0/phase0_feasibility.md`](results/phase0/phase0_feasibility.md) | — |
| `analysis_plan.md` | [`results/phase0/analysis_plan.md`](results/phase0/analysis_plan.md) | v1.0 |

Resultados auxiliares em `results/phase0/`: ajustes alométricos (`model_fits.tsv`), sensibilidades, estados ancestrais, vizinhos filogenéticos de *Melopsittacus*, simulação de poder e figuras.

Fase 0b (`results/phase0b/`): matriz espécie × dados, poder com o conjunto real de espécies, comparação dos dois desfechos, pares irmãos com genoma, validação Python × R e mudanças de regime (PhylogeneticEM).

## Reproduzir

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/00_fetch_external.sh              # clona fontes em commits fixos (config/sources.tsv)
.venv/bin/python scripts/01_build_phase0_tables.py
.venv/bin/python scripts/02_phase0_analysis.py
.venv/bin/python scripts/03_sra_parrot_inventory.py   # ~1 GB de leitura parcial do Parquet do SRA na AWS
.venv/bin/python scripts/04_build_data_matrix.py
.venv/bin/python scripts/05_phase0b_analysis.py
bash scripts/setup_r_envs.sh                          # R via conda-forge (micromamba)
.renv/bin/Rscript scripts/R/validate_phase0.R
.renv45/bin/Rscript scripts/R/regime_shifts.R
```

As análises filogenéticas principais estão em Python (`scripts/lib/phylo.py`) e foram **validadas contra `nlme`, `phylolm` e `phytools`**, com resultados idênticos (`results/phase0b/validation_python_vs_R.json`).

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

- A política de rede bloqueou NCBI (Datasets/E-utilities/FTP), Crossref, Europe PMC, Dryad, Zenodo, Figshare, bioRxiv e CRAN. Os dados de sequenciamento por espécie foram obtidos dos metadados oficiais do SRA no AWS Open Data; R veio do conda-forge. Continuam pendentes: accessions e BUSCO das montagens, a árvore de Smith et al. e a busca bibliográfica sistemática.
- Linhas marcadas com `to_verify` ou `unverified` **não** devem ser usadas como fato.
