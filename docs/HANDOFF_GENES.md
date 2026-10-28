# Handoff: etapa de genes (H8 e H6)

Atualizado em 28 out. 2026. Para quem continuar (Angelo, Gemini/Antigravity ou outro agente). Nenhum gene real foi analisado até esta data.

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
