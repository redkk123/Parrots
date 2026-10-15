# Auditoria das montagens (Fase 1) — relatório parcial

Iniciado em 14 out. 2026, atualizado em 15 out. 2026. **Nenhuma análise de genes foi feita.** Item 1.1 do `TODO.md`.

Estado: inventário NCBI ✅ · identidade mitocondrial 96/96 ✅ · checagem nuclear (Mash) ⏳ · BUSCO 3/96 ⏳ (≈ 1h20 por genoma nesta máquina; ~5 dias para o lote).

## 1. O que é o "lote Hains"

O NCBI tem **316 montagens de Psittaciformes**; **235 são da Iridian Genomes** (202 atuais, 196 espécies, depositadas de 2020 a 2025). O artigo de Hains et al. 2022 cobre só uma parte (94 spp.). Na matriz da Fase 0b, `in_hains2022_iridian` quer dizer, na prática, "tem montagem da Iridian". Das 105 espécies do conjunto de análise, **96** têm montagem Iridian; 9 só têm montagens de outras fontes (VGP, B10K, Sanger etc.), e uma (*Pseudeos fuscata*) não tem montagem no NCBI, só reads.

Tabela: [`assembly_inventory.tsv`](assembly_inventory.tsv) (script `scripts/09_ncbi_assembly_inventory.py`). As accessions pendentes de `data/curated/assemblies_literature.tsv` foram resolvidas com os mesmos metadados.

## 2. Contiguidade: o critério pré-registrado corta mais da metade

Todas as 96 são Illumina + SPAdes, de um único indivíduo. N50 de contig: mediana 14 kb (mín. 256 pb, máx. 66 kb).

| Sinal | Montagens |
| --- | --- |
| N50 de contig < 5 kb | 31 |
| Tamanho fora de 1,0–1,35 Gb (de 0,08 a 2,13 Gb) | 16 |
| GC fora de 40–44,5 % | 7 |
| Scaffold N50 > 1 Mb só com reads curtos (scaffolding guiado por referência) | 3 |

A qualidade piora com a safra: N50 < 5 kb em 4/31 montagens "SPAdes v. 2021", contra 13/17 "SPAdes v. 2024".

O plano (§4) exige **BUSCO ≥ 85 % e N50 de contig ≥ 20 kb**. Pelo N50:

| Corte de N50 (melhor montagem por espécie) | Espécies | Psittaculidae / restante |
| --- | --- | --- |
| ≥ 20 kb (pré-registrado) | **47** de 105 (só 22 das 96 Iridian) | **12** / 35 |
| ≥ 10 kb | 76 | 21 / 55 |
| ≥ 5 kb | 82 | 24 / 58 |
| sem corte | 105 | 34 / 71 |

Com 20 kb, Psittaculidae cai de 34 para 12 espécies e os lóris ficam reduzidos a dois *Trichoglossus*. Isso esvazia as análises "dentro de Psittaculidae" e a comparação periquito × Loriini da H8. A condição de escala da Fase 0 (~60 espécies) não é atendida. **Decisão pendente (do pesquisador):** manter 20 kb ou registrar outro critério antes de olhar genes. O BUSCO, quando terminar, deve pesar nessa decisão.

## 3. BUSCO (compleasm, aves_odb10) — 3 de 96

| Espécie | N50 de contig | Completo | Fragmentado | Ausente |
| --- | --- | --- | --- | --- |
| *Psittacus timneh* | 38,6 kb | 87,7 % | 8,8 % | 3,4 % |
| *Psittacara leucophthalmus* | 19,2 kb | 76,1 % | 16,4 % | 7,2 % |
| *Psittacara erythrogenys* | 0,76 kb | **2,0 %** | 8,4 % | 89,6 % |

Mesmo uma montagem com N50 próximo de 20 kb fica abaixo de 85 %. Montagens com N50 < 1 kb são inutilizáveis para extrair genes.

## 4. Identidade do espécime e contaminação (mitogenoma)

Método: contigs mitocondriais localizados com COX1/CYTB/ND2 de cinco papagaios distantes e comparados a 6.700 sequências de COX1, CYTB e ND2 do GenBank, **excluídas** as depositadas por Hains/Pirro/Iridian (`scripts/11_mito_identity.py`). Todas as 105 espécies têm referência independente.

| Resultado | Montagens |
| --- | --- |
| Só o mitogenoma da espécie declarada | 22 |
| Espécie declarada **+ mitogenoma de outra espécie** | 39 |
| Só mitogenoma de outra espécie ("conflito") | 22 |
| Só cópias nucleares divergentes (NUMTs) | 8 |
| Ambíguo / nada encontrado | 4 / 1 |

**61 das 96 montagens contêm o mitogenoma de outra espécie**, em geral completo (contigs de 17–25 kb). Das 75 espécies "doadoras", **88 % foram sequenciadas pela própria Iridian**, contra 55 % esperados entre as espécies com referência no GenBank. O padrão aponta para **contaminação cruzada entre amostras no laboratório** (bibliotecas no mesmo lote, *index hopping*). Há pares recíprocos (*Anodorhynchus leari* ↔ *Eolophus roseicapilla*; *Poicephalus gulielmi* ↔ *Diopsittaca nobilis*).

Como ler isso:
- O mitogenoma tem centenas a milhares de cópias por célula. Cerca de 1 % de reads contaminantes basta para montar um mitogenoma inteiro e quase não gera contigs nucleares. **Mitogenoma estranho não implica genoma nuclear contaminado**, mas mostra que a contaminação existe e precisa ser medida.
- Um "conflito" pode ser troca de amostra ou só falha em montar o mitogenoma próprio. Casos entre famílias (*Cacatua alba* com *Amazona*, *Guaruba* com *Psittacus*, *Melopsittacus* com *Ara severus*) não são explicáveis por taxonomia. Casos dentro do gênero (*Forpus*, *Agapornis*, *Psittacara*) podem refletir referência mal identificada ou introgressão.
- Limite: registros do GenBank feitos por terceiros a partir dos reads da Iridian não são excluídos pelo filtro de autores.

**Próximo passo, já em execução:** distâncias Mash nucleares (`scripts/16_mash_batch.sh`, `scripts/17_mash_identity.py`) para as 96 montagens, 48 montagens independentes das mesmas espécies e 32 doadoras. Um "conflito" é resolvido como identidade correta se o genoma nuclear fica mais perto do gênero declarado que da doadora.

## 5. Qualidade × longevidade (TODO 1.1, item 4)

PGLS (λ de Pagel) da longevidade relativa sobre cada métrica (`quality_vs_longevity.tsv`). Com as métricas do NCBI, **nenhuma associação**: N50 p = 0,25 (expectativa de vida, n = 96) e 0,20 (longevidade máxima, n = 69). O desvio de tamanho fica no limite só para longevidade máxima (p = 0,08). Repetir com BUSCO e com o conjunto final.

## 6. Desvios e decisões registrados

- BUSCO via compleasm, com aves_odb10 explícito. Registrado no histórico de `results/phase0/analysis_plan.md` (15 out. 2026).
- Critérios do §4 **não** alterados.
- Os genomas não ficam em disco: cada montagem é baixada, analisada e apagada. Ficam versionados os resumos (`busco/*/summary.txt`, `mito/*.json`, tabelas).

## 7. Como reproduzir / retomar

```bash
.venv/bin/python scripts/09_ncbi_assembly_inventory.py   # inventário NCBI
bash scripts/12_audit_batch.sh                            # BUSCO + identidade, lote retomável (dias)
bash scripts/15_identity_batch.sh                         # só identidade mitocondrial (~6 h)
bash scripts/16_mash_batch.sh                             # sketches Mash (~9 h)
.venv/bin/python scripts/17_mash_identity.py              # identidade nuclear
.venv/bin/python scripts/14_audit_summary.py              # tabela consolidada: audit_table.tsv
.venv/bin/python scripts/13_quality_vs_longevity.py
```

Ferramentas em `.bioenv/` (micromamba, bioconda): compleasm 0.2.10, miniprot, hmmer, BLAST+, seqkit, mash, ncbi-datasets-cli. Os lotes param se a máquina suspender. Basta relançá-los, porque eles pulam o que já foi feito.
