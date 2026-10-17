# Handoff para o Antigravity (17 out. 2026)

Leia primeiro: `HANDOFF.md` (visão geral), `TODO.md` (tarefas), `results/phase0/analysis_plan.md` (regras e Histórico) e `results/phase1_audit/audit_report.md`. Dono: Angelo (UNICEPLAC, Brasília), fala português informal. Repo: `redkk123/Parrots`, branch `claude/new-session-lwh7xh`, pasta local `/home/red/parrots`.

## Regras do projeto (não quebrar)
- **Nenhuma análise de genes foi feita.** Listas de genes (módulos M1–M8) e testes estão congelados em `analysis_plan.md` §10. Mudanças só com entrada datada no Histórico, **antes** de olhar genes.
- Não afirmar ineditismo nem "redução em Psittaculidae" como estabelecida: é hipótese sensível à árvore.
- E-mails que o Angelo envia: curtos e naturais, sem cara de IA. Nunca inventar afiliação. Não enviar nada por ele sem pedir.
- **Não gravar o código de acesso do Drury-Reavill no repositório** (está no Gmail do Angelo, e-mail de H. Beaufrère de 17 out., vale 15 dias, até ~1 nov.). Passar por variável de ambiente.

## Estado
| Frente | Estado |
| --- | --- |
| Auditoria dos 96 genomas Iridian | identidade mitocondrial 96/96 e nuclear (Mash) 96/96 prontas; **BUSCO 36/96** (≈ 1 h/genoma, lote `scripts/12_audit_batch.sh` + `scripts/autosave_audit.sh`, que commita/pusha de hora em hora) |
| Excluídas por identidade | 4: "A. nigrigenis" GCA_035041115.1 (é passeriforme), C. haematuropygia, T. haematodus (misturadas), A. jandaya (é Poicephalus). Ver `results/phase1_audit/nuclear_identity_review.tsv` |
| Critério de genomas | revisto em 16 out.: BUSCO ≥ 70 % + gene completo por espécie; o original (≥ 85 % e N50 ≥ 20 kb) é sensibilidade. Pares irmãos pela árvore de Smith 2024 (recalcular quando o BUSCO terminar) |
| Árvore de Smith 2024 | feita (`results/phase1_trees/`): efeito de Psittaculidae −25 %, p = 0,05 (supertree: −31 %, p = 0,003); PGLS com erro de medida validado |
| Câncer | Kapsetaki 2024: periquito 20,8 % × Loriini 1,5 %, mas calopsita doméstica 11,4 % (domesticação é explicação concorrente). Relatório em `results/phase1_cancer/cancer_report.md` |
| Drury-Reavill (UC Davis) | **destravado em 17 out.** (ver abaixo) |

## Pendente, em ordem
1. **Drury-Reavill:** Beaufrère disse para identificar necropsias pelo **número de tecidos** (> 2 ou 3). Testado: `tissues_min=2` → 12.198 casos de papagaios, `=3` → 10.152 (o artigo dele usou 11.737). `scripts/22_drury_reavill_counts.py [tissues_min]` já foi trocado para esse filtro (usa só contagens agregadas da API pública, cache em `data/external/drury_reavill/`), mas **não foi rodado** (o processo foi morto por falta de memória). Rodar com 3 (principal) e 2 (sensibilidade); o arquivo antigo `drury_reavill_necropsy.tsv` (filtro "Post Mortem", não interpretável) foi removido. Depois: taxas por grupo (periquito, calopsita, lóris, amazonas, araras, cacatuas, cinzento), padronizadas por idade; atualizar `cancer_report.md`. Com o código dá para exportar registros (`/api/export?...&code=`), mas as contagens bastam.
2. Deixar o BUSCO terminar (~3 dias). Depois: `python scripts/14_audit_summary.py`, escolher um genoma por espécie, **recalcular pares irmãos** na árvore de Smith (`scripts/21_tree_ancestral_pairs.py`) e fixar a lista final de espécies.
3. **Piloto técnico** (ainda não iniciado): 7–9 genomas VGP/GenomeArk + ~100 genes; extrair com miniprot, alinhar códons (MAFFT+PAL2NAL), medir tempo/memória (≤ 12 GB). Não roda testes de seleção.
4. Depois da lista final: RERconverge contra longevidade relativa, RELAX/aBSREL no tronco de Psittaculidae, dentro e fora de Psittaculidae (`TODO.md` §5).
5. Resposta de Boddy/Kapsetaki (idade e tecido por espécie) ainda não chegou.

## Armadilhas desta máquina
- **11 GB de RAM, 8 núcleos.** O compleasm usa ~4–7 GB; rodar outra coisa pesada junto causa kill por memória (exit 137). Não rodar o BUSCO em paralelo com jobs grandes.
- Ferramentas em `.bioenv/` (micromamba: compleasm, miniprot, hmmer, BLAST+, seqkit, mash, datasets). Python em `.venv/`; R em `.renv45/`. Nenhum está no git.
- Os lotes são **retomáveis**: se o PC reiniciar, relançar `nohup setsid bash scripts/12_audit_batch.sh &` e `scripts/autosave_audit.sh &`.
- Cuidado com `pkill -f`: casa com o próprio shell. Usar `pgrep -f '[x]xx'` com colchete.
- Git local já configurado com o nome/e-mail do Angelo (só neste repo). Commits terminam com `Co-Authored-By`.
- O NCBI Datasets e Zenodo funcionam; o Dryad bloqueia download por script (usar o espelho Zenodo).
