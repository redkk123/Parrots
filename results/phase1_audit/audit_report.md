# Auditoria das montagens (Fase 1) — relatório parcial

Iniciado em 14 out. 2026, atualizado em 28 out. 2026. **Nenhuma análise de genes foi feita.** Item 1.1 do `TODO.md`.

Estado: inventário NCBI ✅ · identidade mitocondrial 96/96 ✅ · identidade nuclear (Mash) 96/96 ✅ · BUSCO Iridian 93/96 ⏳ · BUSCO de 30 montagens de outras fontes ⏳ (na fila, §7). Cerca de 1 h por genoma; salvamento automático no GitHub de hora em hora (`scripts/autosave_audit.sh`).

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

Com 20 kb, Psittaculidae cai de 34 para 12 espécies e os lóris ficam reduzidos a dois *Trichoglossus* (e um deles, *T. haematodus*, foi depois excluído por mistura de genomas; §4). Isso esvazia as análises "dentro de Psittaculidae" e a comparação periquito × Loriini da H8. A condição de escala da Fase 0 (~60 espécies) não é atendida. **Decisão pendente (do pesquisador):** manter 20 kb ou registrar outro critério antes de olhar genes. O BUSCO, quando terminar, deve pesar nessa decisão.

## 3. BUSCO (compleasm, aves_odb10) — 10 de 96 (amostra aleatória do lote)

| Espécie | N50 de contig | Completo | Fragmentado | Ausente |
| --- | --- | --- | --- | --- |
| *Psittacus timneh* | 38,6 kb | **87,7 %** | 8,8 % | 3,4 % |
| *Psittacara leucophthalmus* | 19,2 kb | 76,1 % | 16,4 % | 7,2 % |
| *Pyrrhura leucotis* | 13,9 kb | 74,7 % | 17,1 % | 7,9 % |
| *Agapornis personatus* | 15,9 kb | 73,9 % | 16,8 % | 9,1 % |
| *Loriculus philippensis* | 1,8 kb | 51,1 % | 26,0 % | 22,5 % |
| *Enicognathus leptorhynchus* | 3,4 kb | 43,5 % | 27,3 % | 28,7 % |
| *Poicephalus cryptoxanthus* | 1,2 kb | 42,8 % | 25,0 % | 31,9 % |
| *Eolophus roseicapilla* | 0,33 kb | 2,7 % | 13,5 % | 83,6 % |
| *Psittacara erythrogenys* | 0,76 kb | 2,0 % | 8,4 % | 89,6 % |
| *Coracopsis nigra* | 0,31 kb | 0,9 % | 3,7 % | 95,3 % |

Só 1 de 10 passa no critério de 85 %. Montagens com N50 de 14–19 kb ficam em 74–76 %, com muitos genes **fragmentados** (16–17 %), não ausentes. O corte de BUSCO pré-registrado tende a ser **mais restritivo** que o de N50. Como a ordem do lote é aleatória, estes 10 são uma amostra sem viés do restante.

**Atualização (93 de 96, 28 out.):** mediana de 76 % de BUSCO completo. 55 passam no critério revisto (≥ 70 %); 11 chegariam a 85 % e 10 passam no critério original (85 % e N50 ≥ 20 kb). Nenhuma montagem com N50 de contig < 5 kb passa (máximo 61,5 %); com N50 ≥ 5 kb, a mediana é 80 %. Tabela completa: [`audit_table.tsv`](audit_table.tsv).

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

### Identidade nuclear (Mash)

Distâncias Mash (k = 21, 50 mil hashes) entre 176 montagens: as 96 Iridian, 48 montagens independentes das mesmas espécies e 32 da Iridian das espécies "doadoras" (`scripts/16_mash_batch.sh`, `scripts/17_mash_identity.py`; [`mash_identity.tsv`](mash_identity.tsv)).

Calibração, com montagens de N50 ≥ 5 kb: mesma espécie 0,002–0,009; dentro do gênero, mediana 0,006 (*Agapornis*) e 0,011 (*Cacatua*); entre famílias (*Agapornis* × *Cacatua*) 0,08; a maior distância dentro da ordem (*Nestor* × demais) é 0,10. **Com N50 < 5 kb, a distância é inflada pela fragmentação** (0,01–0,08 para a mesma espécie) e não serve para identidade.

- **Dos 22 "conflitos" mitocondriais**, 9 têm genoma nuclear do gênero declarado, mais perto dele que da doadora (*Guaruba*, *Myiopsitta*, *Melopsittacus*, *Cacatua alba*, *Lorius lory*, *Psittacula derbiana* etc.). Neles o mitogenoma estranho era contaminação, não troca de amostra. 10 não são testáveis por esse critério (doadora congênere ou montagem fragmentada), embora *Amazona aestiva* fique a 0,009 da montagem independente da mesma espécie. Em 2, a doadora não tem montagem. 1 falha (*Aratinga jandaya*, abaixo).
- **Quatro montagens falham** ([`nuclear_identity_review.tsv`](nuclear_identity_review.tsv)):

| Montagem | Diagnóstico | Evidência |
| --- | --- | --- |
| "*Agapornis nigrigenis*" GCA_035041115.1 | **não é papagaio** | ≥ 0,145 de qualquer papagaio (mais que *Nestor* × demais); 0,16 da montagem Sanger da mesma espécie; mitogenoma de 13,7 kb mais próximo de *Garrulax* (Leiothrichidae, Passeriformes), 89,5 % |
| "*Cacatua haematuropygia*" GCA_046039615.1 | mistura com o mesmo genoma estranho | 1,77 Gb; 0,116 das outras *Cacatua* (esperado ~0,01); vizinho mais próximo é a montagem acima |
| "*Trichoglossus haematodus*" GCA_032361345.1 | mistura com o mesmo genoma estranho | 0,128 dos congêneres; vizinho mais próximo é a montagem acima; sem mitogenoma de papagaio |
| "*Aratinga jandaya*" GCA_026284485.1 | **espécime errado** | 0,0039 de *Poicephalus meyeri* (nível de mesma espécie) contra 0,066 do *Aratinga* mais próximo; só mitogenomas de *Poicephalus* |

Para *A. nigrigenis* há montagem cromossômica independente (Sanger, GCA_988293565.1), que deve substituir a da Iridian. As outras três espécies ficam sem genoma utilizável.

A mitocôndria estranha em 61 montagens e esses quatro casos indicam **controle de amostras falho no lote Iridian**. A identidade de cada montagem tem de ser checada antes do uso, como aqui; o metadado do NCBI não basta.

## 5. Qualidade × longevidade (TODO 1.1, item 4)

PGLS (λ de Pagel) da longevidade relativa sobre cada métrica (`quality_vs_longevity.tsv`). Com as métricas do NCBI, **nenhuma associação**: N50 p = 0,25 (expectativa de vida, n = 96) e 0,20 (longevidade máxima, n = 69). O desvio de tamanho fica no limite só para longevidade máxima (p = 0,08). Repetir com BUSCO e com o conjunto final.

## 6. Desvios e decisões registrados

- BUSCO via compleasm, com aves_odb10 explícito. Registrado no histórico de `results/phase0/analysis_plan.md` (15 out. 2026).
- Critérios do §4 revistos em 16 e 28 out. 2026, antes de olhar genes (Histórico de `results/phase0/analysis_plan.md`).
- Os genomas não ficam em disco: cada montagem é baixada, analisada e apagada. Ficam versionados os resumos (`busco/*/summary.txt`, `mito/*.json`, tabelas).

## 7. Um genoma por espécie (provisório, 28 out.)

Regras do plano (Histórico, 16 e 28 out.): BUSCO ≥ 70 %, identidade não reprovada e, entre as montagens aprovadas, preferência pelas de escala cromossômica de outras fontes. Quando o mitogenoma não decide, a identidade é aceita no nível de clado pelo Mash (resíduo em relação à árvore ≤ 0,006, N50 ≥ 5 kb). Isso vale para 12 montagens Iridian, 10 delas aprovadas no BUSCO. Script: `scripts/26_select_genomes.py`; tabela: [`species_genomes.tsv`](species_genomes.tsv).

| Situação (105 espécies) | Loriini | *Melopsittacus* | Outros Psittaculidae | Resto da ordem | Total |
| --- | --- | --- | --- | --- | --- |
| Incluídas | 4 | 0 | 7 | 40 | **51** |
| Aguardando BUSCO (Iridian ou outra fonte) | 0 | 1 | 9 | 13 | 23 |
| Excluídas | 5 | 0 | 7 | 18 | 30 |
| Sem montagem (*Pseudeos fuscata*) | 1 | 0 | 0 | 0 | 1 |

- Das 51 incluídas, 41 têm identidade confirmada (entram na sensibilidade de identidade) e 8 passam no critério original.
- Quatro das incluídas (*Amazona aestiva*, *Ara ararauna*, *Guaruba guarouba*, *Myiopsitta monachus*) podem trocar de genoma: há montagem cromossômica de outra fonte na fila.
- **Fila de outras fontes** ([`audit_queue_other.tsv`](audit_queue_other.tsv), 30 montagens, Psittaculidae primeiro). Entra uma montagem quando acrescenta a espécie, é cromossômica ou tem N50 mais de 2× o da Iridian. Começa sozinha quando o lote Iridian termina (serviço `parrots-busco-other`).
- **Se todas passarem:** até 74 espécies, 21 delas Psittaculidae, incluindo o periquito (VGP), *Lathamus*, *Neophema*, *Neopsephotus*, *Cyclopsitta*, *Psittacula krameri*, *Loriculus galgulus* e três *Agapornis* cromossômicos. Com isso volta a condição de escala da Fase 0 (~60 espécies).
- **Pares irmãos** (árvore de Smith, só espécies incluídas): 17 pares, 10 também irmãos na *supertree* ([`../phase1_trees/sister_pairs_audited.tsv`](../phase1_trees/sister_pairs_audited.tsv)). São provisórios: periquito × Loriini ainda não aparece, porque os genomas dele estão na fila.

## 8. Como reproduzir / retomar

```bash
.venv/bin/python scripts/09_ncbi_assembly_inventory.py   # inventário NCBI
bash scripts/keep_audit_running.sh                        # relança o lote se parou; procura genomas novos no NCBI 1x/semana
bash scripts/12_audit_batch.sh                            # BUSCO + identidade, lote retomável (dias)
bash scripts/15_identity_batch.sh                         # só identidade mitocondrial (~6 h)
bash scripts/16_mash_batch.sh                             # sketches Mash (~9 h)
.venv/bin/python scripts/17_mash_identity.py              # identidade nuclear
.venv/bin/python scripts/14_audit_summary.py              # tabela consolidada: audit_table.tsv
.venv/bin/python scripts/26_select_genomes.py              # um genoma por espécie + fila de outras fontes
bash scripts/12_audit_batch.sh results/phase1_audit/audit_queue_other.tsv   # BUSCO + identidade dessa fila
.venv/bin/python scripts/21_tree_ancestral_pairs.py        # pares irmãos entre as espécies incluídas
.venv/bin/python scripts/13_quality_vs_longevity.py
```

Ferramentas em `.bioenv/` (micromamba, bioconda): compleasm 0.2.10, miniprot, hmmer, BLAST+, seqkit, mash, ncbi-datasets-cli. Os lotes param se a máquina suspender. Basta relançá-los, porque eles pulam o que já foi feito.
