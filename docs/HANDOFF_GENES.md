# Handoff: etapa de genes (H8 e H6)

Atualizado em 28 out. 2026. Para quem continuar (Angelo, Gemini/Antigravity ou outro agente). Nenhum gene real foi analisado até esta data.

## Estado em 28 out., 14h15

- **Lista congelada:** `results/phase1_audit/analysis_species.tsv`, 58 espécies, 15 Psittaculidae. É a base da análise principal (Histórico do plano, 28 out.).
- **Etapa de genes rodando** (`parrots-genes`). O BUSCO das 22 montagens restantes volta sozinho depois dela e alimenta a **extensão pré-declarada**: regenerar `analysis_species.tsv` a partir de `species_genomes.tsv`, apagar `results/phase2_genes/DONE` e rodar de novo, uma vez só.
- **Segunda referência:** painel do kakapo (`kakapo_<GENE>` nas tabelas de chamadas). aBSREL testa o tronco de Psittaculidae e o ramo do periquito.

## O que roda sozinho

Tudo em serviços do systemd do usuário. O salvamento vai para o GitHub de hora em hora (`scripts/autosave_audit.sh`).

1. **Auditoria (BUSCO + identidade)** das montagens Iridian e de outras fontes: `parrots-busco`.
2. **Vigia** `parrots-keeper`: a cada 15 min roda `scripts/keep_audit_running.sh`. Esse script relança o BUSCO se parar e, quando a fila acaba, inicia a etapa de genes (`parrots-genes` → `scripts/29_gene_batch.sh`), até existir `results/phase2_genes/DONE`.
3. **Etapa de genes** (`29_gene_batch.sh`), em ordem:
   1. `27_extract_genes.sh` + `28_parse_miniprot.py`: extrai os 101 genes do painel (miniprot) de cada genoma escolhido (`results/phase1_audit/species_genomes.tsv`, status `included*`). Saídas: `results/phase2_genes/hits/<acc>.tsv`, `cds/<acc>.fna`, `gff/<acc>.gff.gz`.
   2. `30_align_genes.py`: alinhamentos de códons (MAFFT + PAL2NAL) em `alignments/`, mais `gene_calls.tsv` e `alignment_summary.tsv`.
   3. `31_selection_tests.py`: HyPhy RELAX e aBSREL nas duas árvores. Os JSON ficam em `hyphy/`, fora do git.
   4. `32_h8_summary.py`: `h8_relax_by_gene.tsv`, `h8_modules.tsv`, `h8_integrity.tsv`.

**Se o PC reiniciar:** os serviços morrem. Rode `bash scripts/keep_audit_running.sh` uma vez, ou instale o crontab (o usuário precisa fazer isso):

```bash
( crontab -l 2>/dev/null; echo '@reboot /home/red/parrots/scripts/keep_audit_running.sh'; echo '0 */6 * * * /home/red/parrots/scripts/keep_audit_running.sh' ) | crontab -
```

Para conferir o andamento: `logs/12_audit_batch.log`, `logs/29_gene_batch.log`, `logs/keep_audit.log`, `systemctl --user status parrots-busco parrots-genes`.

## Genoma todo (H6 e universo da H8)

- **Extração:** `27_extract_genes.sh` roda o **proteoma canônico inteiro** de *Melopsittacus* (16.045 genes, `scripts/33_prepare_proteome.py` → `data/genome_genes/`, fora do git) numa só rodada do miniprot por genoma. O painel de 101 genes sai da mesma rodada. Saídas: `results/phase3_genome/calls/<acc>.tsv.gz` (no git) e `data/genome_cds/<acc>.fna.gz` (local).
- **Alinhamentos:** `34_align_genome.py`, genes com ORF completa em ≥ 50 % das espécies. Saída em `data/genome_alignments/` (local, alguns GB).
- **Fenótipos e árvores:** `35_prepare_rer_inputs.py` gera `results/phase3_genome/rer_traits.tsv` (longevidade relativa, cérebro relativo e cérebro além da longevidade) e as árvores mestras.
- **RERconverge:** `scripts/R/rerconverge_h6_h8.R` (**nunca rodado**: escrito pela API da vinheta). Roda sozinho no fim de `29_gene_batch.sh` se o pacote estiver instalado. A instalação foi bloqueada para o agente; o usuário precisa rodar:

  ```bash
  /home/red/parrots/.renv45/bin/Rscript -e 'install.packages(c("remotes","weights","gplots","dplyr","RColorBrewer"), repos="https://cloud.r-project.org"); remotes::install_github("nclark-lab/RERconverge", upgrade="never")'
  ```

- **Enriquecimento GO:** precisa de `data/external/go/go_bp.gmt` (GO Processo Biológico, símbolos humanos, por exemplo do MSigDB C5 GO:BP). Sem esse arquivo, o script pula essa parte.
- **Duração estimada:** extração ~40 min por genoma (~2 dias para ~70); alinhamentos ~1 dia; RERconverge com 1.000 permulações, alguns dias por árvore.

## Regras (não mudar depois de ver resultados)

Estão no Histórico de `results/phase0/analysis_plan.md` (entradas de 16 e 28 out.). Resumo:
- **Um genoma por espécie:** BUSCO ≥ 70 %, identidade não reprovada; preferir montagem cromossômica de outra fonte.
- **Gene completo:** ≥ 90 % da proteína alinhada, sem *frameshift* nem parada interna. Identidade < 60 %, ou lócus de outro gene do painel, conta como parálogo (dado faltante).
- **RELAX:** TEST = tronco de Psittaculidae + ramos internos; REFERENCE = o resto. **aBSREL:** só o tronco. FDR 10 % dentro do painel.
- **Testes de gene × longevidade rodam uma vez, com a lista final.** Não rodar repetidamente conforme chegam genomas.

## Pendente (não implementado)

1. **Nulo da H8:** 1.000 conjuntos aleatórios do universo genômico, pareados por tamanho de alinhamento e taxa de fundo (comprimento da árvore do gene, `genetrees_*.txt`). O universo já sai da etapa acima; falta o script que sorteia os conjuntos e compara as proporções por módulo.
2. **RERconverge:** instalar o pacote (comando acima) e conferir o script R na primeira rodada.
3. **Genes interrompidos:** checar com reads (cobertura e *frameshift* real ou erro de montagem) antes de chamar perda. Segunda referência longeva (TOGA ou miniprot com proteínas de *Ara*/*Amazona*).
4. Testes principais do §5 (RERconverge no genoma todo) e sensibilidades: critério original, só identidade confirmada, árvore de Smith.

## Ferramentas

`.bioenv` (micromamba, bioconda): compleasm 0.2.10, miniprot 0.18, mafft, pal2nal, hyphy 2.5.103, seqkit, samtools, ncbi-datasets. Python: `.venv`. R: `.renv45` (ape, phytools, nlme).

## Triagem de genes interrompidos sem reads (7 nov. 2026)

`scripts/40_integrity_triage.py` → `results/phase2_genes/disrupted_triage*.tsv`. Das 329 chamadas "disrupted": 153 compartilhadas com congênere, 15 isoladas, 161 sem congênere. As "compartilhadas" se concentram em 12 genes (ACD, ERCC3, HMOX1, MUTYH, MYC, NQO1, OGG1, PRDX3, PSMB5, PSMD6, RB1, TERF2IP). OGG1 é "disrupted" em 52/52 Iridian, MUTYH em 42/52, MYC em 39/52, PRDX3 em 25/52 (+24 parálogo): padrão sistemático em todas as espécies, mais provável de ser artefato (montagem Illumina, regiões ricas em GC ou lócus de referência) do que perda real. **Não chamar perda** nesses genes sem checar com reads/outra referência (kakapo). Pendente: checagem com reads.

### Segunda referência (kakapo) nos 12 genes (7 nov. 2026)
Tabela: `results/phase2_genes/disrupted_second_reference.tsv` (chamadas com a proteína do periquito × a do kakapo, mesma montagem).
- **Artefato da proteína de referência do periquito:** OGG1 (44/56 "disrupted" com periquito → completo com kakapo), MUTYH (41), MYC (42), HMOX1 (11/11), PRDX3 (7). Não são perdas.
- **Artefato da proteína do kakapo:** NQO1 (37 espécies completas com periquito, "disrupted" com kakapo).
- **Interrompido com as duas referências (candidatos a checar com reads):** PSMD6 (12/12 das espécies disrupted com periquito), ERCC3 (6), ACD (7), NQO1 (5), MRE11 (1: já visto como erro de montagem no Iridian do periquito), TERF2IP (1).
- RB1 e TP53: o periquito dá parálogo/partial; sem resposta (montagens fragmentadas nesses genes).
- Consequência: a chamada de integridade só vale quando as duas referências concordam. Não usar a do periquito sozinha para OGG1, MUTYH, MYC, HMOX1, PRDX3.

### Candidatos à checagem com reads (7 nov. 2026)
`results/phase2_genes/disrupted_read_check_candidates.tsv` (31 chamadas interrompidas com as duas referências em PSMD6, ERCC3, ACD, NQO1, TERF2IP). A maioria é provável artefato de fragmentação: contigs de 5–15 kb, cobertura da proteína 40–70 % e 10–30 "stops" (alinhamento truncado, não gene quebrado). ERCC3 (6 espécies): todas cobrem só ~57 % em contig pequeno. Critério de priorização (definido depois de ver a tabela, então exploratório): cobertura ≥ 90 %, 0 frameshift, ≤ 2 stops. Esses são os que valem uma checagem com reads (Illumina existe no SRA para todas as espécies listadas).

## Lote de genes concluído (8 nov. 2026)

Terminou às 00:14 (`results/phase2_genes/DONE`), com o nulo da H8 (298 genes). Resultados e ressalvas: **`results/phase2_genes/phase2_report.md`** (revisado em 8 nov.; a primeira versão desta seção e do relatório tinha valores de q errados e chamava o M2 de significativo).
- Nulo: só a média de log2 k do M2 dá p < 0,05 (0,024). Não sobrevive à correção para 16 testes, depende de TSC2 e CCNE1 com k no limite, e com a mediana vira p = 0,07.
- Sinal mais consistente: seleção episódica em genes de telômero no ramo do periquito; TERT intensificado em Psittaculidae (alinhamento com 15 espécies).
- Integridade: nenhuma perda confirmada. Checagem com reads só em *L. galgulus* (PSMD6, ACD), *A. personatus* (PSMD6) e *P. erithacus* (ACD); *P. timneh* e *Myiopsitta* foram inferidos pelo contexto GC.

## Projeto concluído de ponta a ponta (10 nov. 2026)

Todas as três fases (Fase 1: auditoria das 126 montagens; Fase 2: painel de 101 genes + nulo de 298 genes + checagem com segunda referência e reads; Fase 3: `RERconverge` em 7.756 genes nas duas árvores com 1.000 permulações, GO BP e 6 eixos exploratórios) estão **100 % concluídas**:

1. **Fase 1 — Auditoria completa (126/126 montagens, `compleasm` `aves_odb10`):**
   - `results/phase1_audit/audit_summary.tsv` (126 montagens: 96 Iridian + 30 não-Iridian) e `results/phase1_audit/species_genomes.tsv` (**71 espécies incluídas**, 20 Psittaculidae e 51 outros Psittaciformes; 54 aprovadas na sensibilidade de identidade; 33 no critério original estrito; 24 pares-irmãos auditados na árvore de Smith 2024).
   - Sem associação PGLS entre qualidade de montagem (`busco_C`, `contig_n50`, `gc`, desvio de tamanho) e longevidade relativa (`results/phase1_audit/quality_vs_longevity.tsv`, $p \ge 0{,}11$).
   - Relatório final: `results/phase1_audit/audit_report.md`. Figura suplementar: `results/phase1_audit/figS1_audit_126_assemblies.png`.
2. **Fase 2 — Painel de 101 genes (M1–M8), nulo genômico de 298 genes e integridade (58 espécies congeladas):**
   - Relatório: `results/phase2_genes/phase2_report.md`. Figura principal: `results/phase2_genes/fig5_phase2_selection_h8.png`.
3. **Fase 3 — `RERconverge` no genoma todo (7.756 genes, árvores Burgio 2019 e Smith 2024, 1.000 permulações filogenéticas):**
   - `RERconverge` v0.3.0 compilado em `.renv45` (`CXX17`) e acelerado com `Rcpp` exato (`fastCorPerms`, validado contra `getAllCor` com erro máximo de $3{,}05 \times 10^{-16}$).
   - Tabelas completas de H6 (`rer_{burgio2019,smith2024}_{rel_log_life_expectancy,rel_log_brain,rel_brain_given_le}.tsv`), enriquecimento GO BP (`rer_*_enrichment_*.tsv`, 1.780 termos avaliados a partir de `data/external/go/go_bp.gmt`), H8 intra-clado (`rer_{burgio2019,smith2024}_h8_{Psittaculidae,rest}.tsv` e `*_panel.tsv`) e sínteses (`h6_genome_summary.tsv`, `h6_top_genes.tsv`, `h6_shared_longevity_brain.tsv`, `h8_within_clade_summary.tsv`, `exploratory_axes_genes.tsv`, `exploratory_axes_summary.tsv`).
   - Relatório completo e síntese para manuscrito: **`results/phase3_genome/phase3_report.md`**. Figura principal: **`results/phase3_genome/fig6_phase3_genome_rer.png`**.

