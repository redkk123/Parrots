# Passagem de bastão: o que é este projeto e onde paramos

Para retomar numa nova sessão (ou no computador local). Leia este arquivo primeiro, depois `TODO.md`. O protocolo completo está em `docs/protocol.md`.

## 1. O objetivo, em uma frase

Entender **como a longevidade evoluiu nos papagaios (Psittaciformes)** e se ela se associa a **mudanças moleculares em genes de manutenção do corpo** (reparo de DNA, controle de câncer, resposta ao dano oxidativo, telômeros), separando isso de ancestralidade compartilhada, tamanho do corpo e qualidade dos dados.

É uma investigação **comparativa e computacional** (genomas e árvores filogenéticas), sem laboratório.

## 2. As três perguntas integradas

1. **Câncer:** por que o periquito-australiano tem mais tumores que outros papagaios (20,6 % vs. 3,6 % numa série clínica de 2002)? Quem são os parentes certos para comparar?
2. **Neurodegeneração:** a raridade de relatos de degeneração em papagaios idosos significa proteção, subdiagnóstico ou falta de estudos?
3. **Radicais livres:** como eles mantêm função com metabolismo alto? Produzem menos radicais, neutralizam mais, reparam melhor ou toleram o dano?

## 3. Sua hipótese (H8)

Papagaios longevos teriam **travas moleculares modulares** (checkpoints, reparo, controle de proliferação, apoptose, resposta ao estresse) que impedem o dano de virar câncer ou disfunção. O periquito teria **só parte** delas.

Foi revisada em 26 set. e dividida em três partes que não se confundem:
- **H8a (proteção modular):** o periquito repara bem o dano oxidativo, mas é menos eficaz em frear células danificadas.
- **H8b (efeito duplo):** parte da resistência vem de células danificadas **sobreviverem mais**, o que também favorece tumores.
- **H8c (história):** as diferenças surgiram no ramo que leva a Psittaculidae (o grupo do periquito).

Regras fixadas **antes** de olhar genes: listas de genes (módulos M1–M8) e testes estão congelados em `results/phase0/analysis_plan.md`, seção 10. Mudar exige entrada datada. Genes intactos enfraquecem só a versão de perda gênica; a versão funcional precisa de dados regulatórios e celulares.

## 4. O que já foi feito (Fases 0 e 0b)

| Resultado | Onde ver |
| --- | --- |
| 218 espécies com expectativa de vida (cativeiro) e massa; 173 com longevidade máxima | `results/phase0/phenotypes.tsv`, `results/phase0b/species_data_matrix.tsv` |
| Alometria: λ de Pagel = 0,65, inclinação 0,35 | `results/phase0/model_fits.tsv` |
| **105 espécies com genoma utilizável** (≥ 20× ou long-read) e fenótipo; 96 vêm do lote de Hains et al. 2022 | `results/phase0b/phase0b_report.md` |
| A medida de longevidade importa: expectativa de vida × recorde máximo concordam pouco (ρ = 0,39) | idem |
| 15 espécies com recorde menor que a média (erro em alguma fonte) | coluna `qc_outcomes_inconsistent` |
| História: ancestral relativamente longevo e **redução em Psittaculidae** (−31 %); **não** há vários aumentos independentes. *Sensível à árvore:* com a árvore de Smith 2024, −25 % e p = 0,05 | `phase0b_report.md` §7; `results/phase1_trees/tree_sensitivity_report.md` |
| Python validado contra R (`nlme`, `phylolm`, `phytools`): resultados idênticos | `results/phase0b/validation_python_vs_R.json` |
| Decisão: **cenário B** (associação contínua, sem afirmar múltiplas origens) | `results/phase0/phase0_feasibility.md` |

**Cuidados já conhecidos:** a conclusão sobre o kea depende do desfecho; o achado de Psittaculidae é **um único evento** (confundido com qualquer outra coisa que mudou ali); a cobertura genômica é estimada, não medida; vários artigos foram lidos só pelo resumo.

## 5. O que fazer agora (resumo; detalhes em `TODO.md`)

1. **Auditar os genomas de Hains et al.** (BUSCO, identidade da espécie pelo mitogenoma, contaminação). É o passo mais importante: 96 das 105 espécies dependem dele.
2. Piloto técnico com os 7–9 genomas VGP (GenomeArk) e ~100 genes.
3. Teste principal: RERconverge contra o fenótipo contínuo e RELAX no tronco de Psittaculidae, repetido **dentro** de Psittaculidae e **dentro** do restante.
4. Eixo do câncer: periquito × Loriini com denominadores de necropsia.
5. Experimento celular H8a × H8b (precisa de laboratório parceiro).

## 6. Como rodar no computador local

```bash
git clone https://github.com/redkk123/Parrots.git
cd Parrots
git checkout claude/new-session-lwh7xh

python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/00_fetch_external.sh                       # baixa fontes em commits fixos
.venv/bin/python scripts/01_build_phase0_tables.py
.venv/bin/python scripts/02_phase0_analysis.py
.venv/bin/python scripts/03_sra_parrot_inventory.py     # ~1 GB de leitura parcial do SRA (AWS)
.venv/bin/python scripts/04_build_data_matrix.py
.venv/bin/python scripts/05_phase0b_analysis.py
bash scripts/setup_r_envs.sh                            # R via conda-forge (micromamba)
.renv/bin/Rscript scripts/R/validate_phase0.R
.renv45/bin/Rscript scripts/R/regime_shifts.R
.venv/bin/python scripts/06_regime_followup.py
```

Os resultados já estão versionados em `results/`, então **não é preciso rodar tudo** para continuar: o que está no repositório basta para começar.

**Diferença importante do ambiente local:** na sessão em nuvem, NCBI, Dryad, Zenodo, Crossref e bioRxiv estavam bloqueados. No seu computador isso provavelmente não vale. Por isso as tarefas marcadas 🔴 no `TODO.md` (accessions e BUSCO no NCBI, árvore de Smith et al. no Dryad, texto completo dos artigos, dados de Kapsetaki) devem ser **mais fáceis localmente**.

Para espaço em disco: reads brutos de 96 genomas short-read são grandes (centenas de GB no total). Comece só pelas montagens do GenomeArk e por um subconjunto.

## 7. Mapa dos arquivos

```
docs/protocol.md            protocolo completo (inclui H1–H8)
TODO.md                     lista de tarefas por prioridade
HANDOFF.md                  este arquivo
config/sources.tsv          fontes externas e versões fixadas
data/curated/               tabelas curadas à mão (literatura, câncer, neuro, oxidativo, nomes)
results/phase0/             Fase 0: fenótipos, árvore, análises, plano de análise
results/phase0b/            Fase 0b: matriz espécie × dados, regimes, relatório
scripts/                    pipeline numerado (00–06), lib/ (Python) e R/
logs/                       logs e checksums das fontes
```

## 8. Para colar na nova sessão

> Estou retomando um projeto de genômica comparativa da longevidade em papagaios. Leia `HANDOFF.md`, `TODO.md` e `results/phase0b/phase0b_report.md`. A hipótese H8 e os testes já estão fixados em `docs/protocol.md` e `results/phase0/analysis_plan.md` (seção 10); não os altere sem uma entrada datada. Começamos pela auditoria dos genomas de Hains et al. (BUSCO, identidade da espécie, contaminação), antes de qualquer análise de genes. Aqui no meu computador a rede não deve estar bloqueada.
