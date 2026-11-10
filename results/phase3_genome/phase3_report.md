# Fase 3 — Correlações Genômicas de Taxa Evolutiva (RERconverge: H6, H8 Intra-Clado e Eixos Exploratórios)

**Data de conclusão:** 10 de novembro de 2026  
**Amostra congelada:** 58 espécies (`results/phase1_audit/analysis_species.tsv`; 15 Psittaculidae e 43 outros Psittaciformes)  
**Universo genômico elegível:** 7.756 genes codificantes com ORF completa em $\ge 50\%$ das espécies (`results/phase3_genome/genome_alignment_summary.tsv`)  
**Árvores avaliadas:** Burgio et al. 2019 (principal, `master_burgio2019.nwk`) e Smith et al. 2024 (sensibilidade, `master_smith2024.nwk`)  
**Permulações filogenéticas:** 1.000 simulações brownianas condicionadas à topologia (`getPermsContinuous` + `fastCorPerms`, semente `20261015`)  
**Script reprodutível:** `scripts/R/rerconverge_h6_h8.R` e `scripts/41_phase3_synthesis.py`  
**Figura síntese:** `results/phase3_genome/fig6_phase3_genome_rer.png`

---

## 1. Visão Geral e Desenho Analítico (H6 e H8 Intra-Clado)

Conforme pré-registrado no plano de análise (`results/phase0/analysis_plan.md`, §5 e entradas de 16 out., 28 out. e 10 nov. 2026), estimamos comprimentos de ramo em substituições de aminoácidos sob o modelo empírico **LG** (`phangorn::optim.pml`) para todos os 7.756 alinhamentos elegíveis nas duas topologias filogenéticas (`genetrees_burgio2019.txt` e `genetrees_smith2024.txt`). Em seguida, extraímos as taxas evolutivas relativas por ramo (`RERconverge::getAllResiduals`, transformação raiz quadrada, ponderada e normalizada) e testamos a associação contínua (`correlateWithContinuousPhenotype` com *winsorization* `winsorizeRER = 3, winsorizeTrait = 3`) contra três fenótipos residuais controlados por massa corporal via PGLS (`results/phase3_genome/rer_traits.tsv`):

1. **`rel_log_life_expectancy`**: longevidade relativa (resíduo PGLS de $\log_{10}(\text{expectativa de vida})$ sobre $\log_{10}(\text{massa corporal})$, $n = 58$ espécies, mediana de 97 ramos por gene).
2. **`rel_log_brain`**: encefalização relativa (resíduo PGLS de $\log_{10}(\text{massa cerebral})$ sobre $\log_{10}(\text{massa corporal})$, $n = 50$ espécies com dado endocraniano, mediana de 85 ramos por gene).
3. **`rel_brain_given_le`**: tamanho cerebral relativo além da longevidade (resíduo PGLS de `rel_log_brain` sobre `rel_log_life_expectancy`, $n = 50$ espécies, mediana de 85 ramos por gene).

Para a **H8 intra-clado**, testamos `rel_log_life_expectancy` separadamente dentro de **Psittaculidae** ($n = 15$ espécies, 7.045 genes com cobertura suficiente, mediana de 23 ramos) e no **restante de Psittaciformes** (`rest`, $n = 43$ espécies, 7.756 genes, mediana de 73 ramos), extraindo o subconjunto do painel pré-registrado (`*_panel.tsv`) com correção de Benjamini-Hochberg restrita ao painel (`panelQ` sobre `permP` e `panel_padj_analytic` sobre `P`).

### 1.1 Resumo Quantitativo Global (`h6_genome_summary.tsv`)

| Árvore | Fenótipo / Recorte | Genes testados | Ramos (mediana) | $q_{\text{analítico}} < 0{,}10$ (`p.adj`) | $q_{\text{perm}} < 0{,}10$ (`permQ`) | $P < 0{,}01$ | $\text{permP} < 0{,}01$ | $P < 0{,}01$ & $\text{permP} < 0{,}05$ | Top gene ($\text{Rho} > 0$) | Top gene ($\text{Rho} < 0$) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Burgio 2019** | `rel_log_life_expectancy` | 7.756 | 97 | 252 | 0 (`min 0,358`) | 322 | 201 | **287** | `YWHAQ` ($+0{,}567$) | `LOC101871916` ($-0{,}525$) |
| **Burgio 2019** | `rel_log_brain` | 7.756 | 85 | 190 | 0 (`min 0,794`) | 286 | 86 | **207** | `MAD2L2` ($+0{,}534$) | `TOMM20` ($-0{,}498$) |
| **Burgio 2019** | `rel_brain_given_le` | 7.756 | 85 | 182 | 0 (`min 0,803`) | 278 | 84 | **194** | `MAD2L2` ($+0{,}534$) | `TOMM20` ($-0{,}486$) |
| **Burgio 2019** | `h8_Psittaculidae` | 7.045 | 23 | 0 | 0 | 80 | 66 | **62** | `USP37` ($+0{,}766$) | `UQCRFS1` ($-0{,}717$) |
| **Burgio 2019** | `h8_rest` | 7.756 | 73 | 293 | 0 (`min 0,257`) | 365 | 259 | **324** | `SKIC8` ($+0{,}636$) | `HNRNPD` ($-0{,}461$) |
| **Smith 2024** | `rel_log_life_expectancy` | 7.756 | 97 | 277 | 0 (`min 0,683`) | 324 | 116 | **241** | `UBL3` ($+0{,}560$) | `NR2F2` ($-0{,}526$) |
| **Smith 2024** | `rel_log_brain` | 7.756 | 85 | 213 | 0 (`min 0,722`) | 295 | 101 | **215** | `MAD2L2` ($+0{,}566$) | `SSU72` ($-0{,}524$) |
| **Smith 2024** | `rel_brain_given_le` | 7.756 | 85 | 205 | 0 (`min 0,750`) | 289 | 98 | **203** | `MAD2L2` ($+0{,}565$) | `SSU72` ($-0{,}518$) |
| **Smith 2024** | `h8_Psittaculidae` | 7.045 | 23 | 0 | 0 | 80 | 64 | **60** | `USP37` ($+0{,}782$) | `DECR2` ($-0{,}753$) |
| **Smith 2024** | `h8_rest` | 7.756 | 73 | 319 | 0 (`min 0,430`) | 389 | 200 | **293** | `SKIC8` ($+0{,}690$) | `HNRNPD` ($-0{,}493$) |

### 1.2 Por que `p.adj` analítico encontra 182–277 genes e `permQ` encontra 0?

Há duas razões complementares e importantes para interpretar essa diferença com rigor:
1. **Estrutura filogenética do fenótipo:** As permulações simulam evolução browniana na árvore de 58 espécies e descontam a pseudorreplicação de clados inteiros (como a divergência basal entre Psittaculidae e Arini/Androglossini). Por isso, `permP` é mais conservador que o $P$ paramétrico.
2. **Piso de resolução de 1.000 permulações em 7.756 genes:** No `RERconverge::permpvalcor`, o $p$-valor empírico bilateral tem piso mínimo $(0 + 1)/(500 + 1) = 1/501 \approx 0{,}001996$. Quando se aplica Benjamini-Hochberg sobre $M = 7.756$ genes, para que o grupo de genes empatados no piso $\text{permP} = 0{,}001996$ alcance $q < 0{,}10$, seriam necessários pelo menos $k \ge 7.756 \times 0{,}001996 / 0{,}10 \approx 155$ genes simultaneamente no piso absoluto ($0$ permulações mais extremas em 500). Na prática, 287 genes (longevidade, Burgio) e 207 genes (cérebro, Burgio) atingem simultaneamente $P < 0{,}01$ e $\text{permP} < 0{,}05$, mas nenhum gene isolado deve ser tratado como "confirmado em FDR genômico de permulação"; em vez disso, a evidência mais sólida vem da **convergência entre longevidade e cérebro nas duas árvores**.

---

## 2. Acoplamento Genômico entre Longevidade e Tamanho Cerebral (H6)

O achado estrutural mais marcante da Fase 3 é o **acoplamento em escala genômica** entre a seleção associada à longevidade relativa e a seleção associada ao tamanho cerebral relativo em Psittaciformes:

- **Correlação entre $\text{Rho}(\text{longevidade})$ e $\text{Rho}(\text{cérebro})$ nos 7.756 genes:**
  - **Burgio 2019:** $r = 0{,}666$ ($p < 10^{-300}$)
  - **Smith 2024:** $r = 0{,}667$ ($p < 10^{-300}$)
- **Persistência no fenótipo residual `rel_brain_given_le` ("cérebro além da longevidade"):**
  - Como a regressão PGLS entre `rel_log_brain` e `rel_log_life_expectancy` na árvore completa tem inclinação modesta ($\beta = 0{,}040 \pm 0{,}032$, $\lambda = 0{,}90$), o resíduo `rel_brain_given_le` preserva quase integralmente a variação de encefalização ($\text{cor}(\text{Rho}_{\text{brain}}, \text{Rho}_{\text{brain|le}}) = 0{,}999$) e mantém forte correlação com a taxa evolutiva associada à longevidade ($r = 0{,}639$ em ambas as árvores).
- **Genes conjuntamente significativos ($P < 0{,}05$ e $\text{permP} < 0{,}05$ com mesmo sinal em longevidade e cérebro; `h6_shared_longevity_brain.tsv`):**
  - **124 genes** na árvore de Burgio 2019 e **123 genes** na árvore de Smith 2024.

### 2.1 Principais Genes Compartilhados entre Longevidade, Cérebro e Supressão Tumoral / Proteostase

| Gene | Função Biológica Principal | $\text{Rho}_{\text{LE}}$ (Burgio / Smith) | $\text{permP}_{\text{LE}}$ (B / S) | $\text{Rho}_{\text{Brain}}$ (Burgio / Smith) | $\text{permP}_{\text{Brain}}$ (B / S) | $\text{Rho}_{\text{Brain\|LE}}$ (B / S) | Direção em espécies longevas e encefalizadas |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **`PTEN`** (painel M2) | Fosfatase supressora tumoral mestre; freio da via PI3K/AKT/mTOR e do crescimento somático/neural | $-0{,}374$ / $-0{,}389$ | $0{,}008$ / $0{,}018$ | $-0{,}424$ / $-0{,}467$ | $0{,}010$ / $0{,}004$ | $-0{,}418$ / $-0{,}462$ | **Forte conservação purificadora ($\text{Rho} < 0$)** |
| **`CTNNB1`** | $\beta$-catenina; eixo canônico Wnt, neurogênese cortical/palial e oncogênese | $-0{,}439$ / $-0{,}437$ | $0{,}004$ / $0{,}012$ | $-0{,}459$ / $-0{,}459$ | $0{,}012$ / $0{,}008$ | $-0{,}447$ / $-0{,}447$ | **Forte conservação purificadora ($\text{Rho} < 0$)** |
| **`CTNNBIP1`** | ICAT; inibidor direto da interação $\beta$-catenina/TCF | $-0{,}362$ / $-0{,}367$ | $0{,}020$ / $0{,}036$ | $-0{,}443$ / $-0{,}444$ | $0{,}004$ / $0{,}008$ | $-0{,}443$ / $-0{,}444$ | **Conservação purificadora ($\text{Rho} < 0$)** |
| **`PSMA4`** (painel M8) | Subunidade $\alpha 4$ do núcleo 20S do proteassoma | $-0{,}397$ / $-0{,}406$ | $0{,}014$ / $0{,}024$ | $-0{,}422$ / $-0{,}434$ | $0{,}014$ / $0{,}008$ | $-0{,}412$ / $-0{,}423$ | **Conservação purificadora ($\text{Rho} < 0$)** |
| **`SSU72`** | Fosfatase do CTD da RNA Pol II e coesão de cromátides irmãs | $-0{,}420$ / $-0{,}463$ | $0{,}006$ / $0{,}002$ | $-0{,}472$ / $-0{,}524$ | $0{,}002$ / $0{,}002$ | $-0{,}466$ / $-0{,}518$ | **Conservação purificadora ($\text{Rho} < 0$)** |
| **`TOMM20`** | Receptor da translocase da membrana externa mitocondrial | $-0{,}325$ / $-0{,}341$ | $0{,}006$ / $0{,}010$ | $-0{,}498$ / $-0{,}509$ | $0{,}002$ / $0{,}002$ | $-0{,}486$ / $-0{,}498$ | **Conservação purificadora ($\text{Rho} < 0$)** |
| **`ANAPC11`** | Subunidade catalítica RING do complexo promotor da anáfase (APC/C) | $-0{,}386$ / $-0{,}402$ | $0{,}014$ / $0{,}014$ | $-0{,}453$ / $-0{,}475$ | $0{,}016$ / $0{,}008$ | $-0{,}447$ / $-0{,}469$ | **Conservação purificadora ($\text{Rho} < 0$)** |
| **`MAD2L2`** (`REV7`) | Subunidade do complexo *shieldin* (reparo de quebras duplas de DNA) e checkpoint do fuso mitótico | $+0{,}418$ / $+0{,}436$ | $0{,}012$ / $0{,}010$ | $+0{,}534$ / $+0{,}566$ | $0{,}002$ / $0{,}002$ | $+0{,}534$ / $+0{,}565$ | **Aceleração evolutiva ($\text{Rho} > 0$, #1 em cérebro)** |
| **`HSPA8`** (`HSC70`) | Chaperona constitutiva neuronal mestre e mediadora da autofagia mediada por chaperonas (CMA) | $+0{,}404$ / $+0{,}399$ | $0{,}008$ / $0{,}012$ | $+0{,}354$ / $+0{,}356$ | $0{,}030$ / $0{,}044$ | $+0{,}349$ / $+0{,}350$ | **Aceleração evolutiva ($\text{Rho} > 0$)** |
| **`RTN1`** | Reticulon-1; morfologia do retículo endoplasmático neuronal e exocitose sináptica | $+0{,}516$ / $+0{,}534$ | $0{,}002$ / $0{,}002$ | $+0{,}440$ / $+0{,}459$ | $0{,}006$ / $0{,}010$ | $+0{,}421$ / $+0{,}440$ | **Aceleração evolutiva ($\text{Rho} > 0$)** |
| **`SOX14`** | Fator de transcrição determinante de interneurônios GABAérgicos e circuitos motores/vocais | $+0{,}476$ / $+0{,}467$ | $0{,}002$ / $0{,}014$ | $+0{,}500$ / $+0{,}488$ | $0{,}002$ / $0{,}004$ | $+0{,}496$ / $+0{,}482$ | **Aceleração evolutiva ($\text{Rho} > 0$)** |
| **`NCALD`** | Neurocalcina delta; sensor neuronal de $\text{Ca}^{2+}$ e endocitose sináptica | $+0{,}417$ / $+0{,}413$ | $0{,}014$ / $0{,}032$ | $+0{,}493$ / $+0{,}501$ | $0{,}002$ / $0{,}006$ | $+0{,}486$ / $+0{,}494$ | **Aceleração evolutiva ($\text{Rho} > 0$)** |
| **`CBLN4`** | Cerebelina-4; organizador trans-sináptico de sinapses inibitórias no pálio/córtex | $+0{,}456$ / $+0{,}476$ | $0{,}004$ / $0{,}004$ | $+0{,}401$ / $+0{,}401$ | $0{,}012$ / $0{,}006$ | $+0{,}388$ / $+0{,}388$ | **Aceleração evolutiva ($\text{Rho} > 0$)** |

**Interpretação biológica:** `PTEN` e `CTNNB1`/`CTNNBIP1` são exemplos clássicos de pleiotropia entre **controle proliferativo (anti-câncer)** e **arquitetura cerebral**. Mutações ou perda de restrição em `PTEN` causam macrocefalia desorganizada e tumorigênese; em psitacídeos longevos e de cérebro grande (como *Ara*, *Cacatua*, *Amazona* e *Strigops*), `PTEN` e o eixo Wnt/$\beta$-catenina sofrem purificação intensificada ($\text{Rho} \approx -0{,}37$ a $-0{,}47$), enquanto reguladores sinápticos finos (`RTN1`, `SOX14`, `NCALD`, `CBLN4`), a chaperona de autofagia neuronal `HSPA8` e a subunidade de reparo de DNA/fuso `MAD2L2` apresentam aceleração evolutiva coordenada.

---

## 3. Enriquecimento Funcional GO Processo Biológico e as 4 Categorias Pré-Registradas da H6

Testamos 1.780 termos GO Biological Process ($\ge 10$ genes presentes entre os 7.756 alinhamentos; `rer_*_enrichment_*.tsv`) via `RERconverge::fastWilcoxGMTAll`. Quando corrigidos para os 1.780 termos simultâneos, nenhum termo isolado atinge $q < 0{,}10$ (`min p.adj = 0,187` para longevidade em Burgio; `0,285` para cérebro em Burgio), porque termos amplos contêm tanto genes sob forte restrição purificadora ($\text{Rho} < 0$) quanto genes acelerados ($\text{Rho} > 0$).

Entretanto, ao examinar as **4 categorias pré-registradas da H6** (`results/phase0/analysis_plan.md`, §5: reparo de DNA/telômeros, ciclo celular/senescência, autofagia/proteostase, mitocôndria/estresse oxidativo) e termos sinápticos, observamos consistência direcional entre as duas árvores em nível nominal ($p < 0{,}05$):

1. **Reparo de DNA e Manutenção Telomérica (Restrição purificadora, $\text{stat} < 0$):**
   - `Regulation Of DNA Recombination (GO:0000018)` ($n = 18$ genes): $\text{stat} = -0{,}245, p = 3{,}2\times 10^{-4}$ (`rel_log_brain` Burgio; $\text{stat} = -0{,}235, p = 5{,}7\times 10^{-4}$ Smith) — **2º termo mais extremo do genoma** para tamanho cerebral.
   - `DNA Damage Response (GO:0006974)` ($n = 149$ genes): $\text{stat} = -0{,}055, p = 0{,}021$ (`rel_log_life_expectancy` Burgio) e $\text{stat} = -0{,}064, p = 0{,}0076$ (`rel_log_brain` Burgio; $p = 0{,}018$ Smith).
   - `Regulation Of Telomere Maintenance (GO:0032204)` ($n = 12$ genes): $\text{stat} = -0{,}203, p = 0{,}015$ (`rel_log_life_expectancy` Burgio; $p = 0{,}013$ Smith) e $\text{stat} = -0{,}208, p = 0{,}013$ (`rel_log_brain` Burgio; $p = 0{,}017$ Smith).
2. **Ciclo Celular e Checkpoints Mitóticos (Restrição purificadora, $\text{stat} < 0$):**
   - `Mitotic G2/M Transition Checkpoint (GO:0044818)` ($n = 25$ genes): $\text{stat} = -0{,}176, p = 0{,}0023$ (`rel_log_brain` Burgio).
   - `Mitotic Cell Cycle Phase Transition (GO:0044772)` ($n = 51$ genes): $\text{stat} = -0{,}111, p = 0{,}0065$ (`rel_log_brain` Burgio; $p = 0{,}0076$ Smith).
3. **Autofagia, Proteostase e Dobramento Proteico (Aceleração evolutiva, $\text{stat} > 0$):**
   - `Negative Regulation Of Protein-Containing Complex Assembly (GO:0031333)` ($n = 26$ genes): $\text{stat} = +0{,}183, p = 0{,}0012$ (`rel_log_life_expectancy` Smith) e $\text{stat} = +0{,}193, p = 6{,}9\times 10^{-4}$ (`rel_log_brain` Burgio).
   - `Regulation Of Proteasomal Ubiquitin-Dependent Protein Catabolic Process (GO:0032434)` ($n = 44$ genes): $\text{stat} = +0{,}098, p = 0{,}025$ (Burgio) e $\text{stat} = +0{,}101, p = 0{,}021$ (Smith).
   - `Macroautophagy (GO:0016236)` ($n = 39$ genes): $\text{stat} = +0{,}099, p = 0{,}032$ (Burgio) e $\text{stat} = +0{,}101, p = 0{,}030$ (Smith).
   - `'De Novo' Post-Translational Protein Folding (GO:0051084)` ($n = 15$ genes): $\text{stat} = +0{,}152, p = 0{,}041$ (`rel_log_life_expectancy` Burgio) e $\text{stat} = +0{,}170, p = 0{,}022$ (`rel_log_brain` Smith).
4. **Tradução e Expressão Mitocondrial (Restrição purificadora, $\text{stat} < 0$):**
   - `Mitochondrial Translation (GO:0032543)` ($n = 57$ genes): $\text{stat} = -0{,}118, p = 0{,}0021$ (`rel_log_life_expectancy` Burgio; $p = 0{,}0077$ Smith).
   - `Mitochondrial Gene Expression (GO:0140053)` ($n = 62$ genes): $\text{stat} = -0{,}105, p = 0{,}0044$ (`rel_log_life_expectancy` Burgio; $p = 0{,}013$ Smith).
5. **Montagem e Reciclagem Sináptica (Aceleração evolutiva, $\text{stat} > 0$):**
   - `Synaptic Vesicle Recycling (GO:0036465)` ($n = 16$ genes): $\text{stat} = +0{,}215, p = 0{,}0030$ (`rel_log_brain` Burgio; $p = 0{,}0051$ Smith).
   - `Positive Regulation Of Synapse Assembly (GO:0051965)` ($n = 12$ genes): $\text{stat} = +0{,}228, p = 0{,}0063$ (`rel_log_life_expectancy` Burgio; $p = 0{,}0086$ Smith).

---

## 4. Os 6 Eixos Funcionais Exploratórios (`exploratory_axes_summary.tsv` e `exploratory_axes_genes.tsv`)

A pedido do usuário e registrado no Histórico de `results/phase0/analysis_plan.md` (10 nov. 2026) antes de rodar a síntese, testamos se seis conjuntos funcionais curados apresentavam deslocamento sistemático de `Rho` em relação ao fundo genômico de 7.756 genes (teste bilateral de Mann-Whitney):

| Eixo Exploratório | Genes presentes (de solicitados) | Mediana $\text{Rho}_{\text{LE}}$ (B / S) | $p_{\text{MW}}$ Longevidade (B / S) | Mediana $\text{Rho}_{\text{Brain}}$ (B / S) | $p_{\text{MW}}$ Cérebro (B / S) | Genes individuais com sinal nominal ($P < 0{,}05$) |
| --- | --- | --- | --- | --- | --- | --- |
| **1. Neuroproteção e Proteostase** | 22 / 32 | $+0{,}039$ / $+0{,}044$ | **$0{,}046$ / $0{,}031$** | $+0{,}068$ / $+0{,}080$ | **$0{,}027$ / $0{,}014$** | `HSPA8` ($\text{Rho} = +0{,}40$, $\text{permP} = 0{,}008\text{--}0{,}012$), `ULK1` ($\text{Rho} = +0{,}27$), `HSP90AB1` ($\text{Rho} = +0{,}25$), `SOD1` ($\text{Rho} = +0{,}22$) |
| **2. Aprendizado Vocal, Linguagem e Plasticidade** | 13 / 25 | $+0{,}009$ / $+0{,}001$ | $0{,}611$ / $0{,}562$ | $-0{,}032$ / $-0{,}018$ | $0{,}193$ / $0{,}351$ | `GAP43` ($\text{Rho}_{\text{LE}} = -0{,}24, P = 0{,}033$); `FOXP2`, `CNTNAP2` e `EGR1` fortemente conservados sem variação de taxa com longevidade |
| **3. Antioxidantes e Equilíbrio Redox** | 22 / 25 | $0{,}000$ / $+0{,}002$ | $0{,}865$ / $0{,}941$ | $+0{,}001$ / $-0{,}002$ | $0{,}809$ / $0{,}780$ | `TXNRD2` ($\text{Rho}_{\text{Brain}} = +0{,}23, \text{permP} = 0{,}026$), `SOD1` ($\text{Rho} = +0{,}22$) |
| **4. Longevidade Extrema e Sensores Nutricionais** | 9 / 24 | $+0{,}039$ / $+0{,}034$ | $0{,}694$ / $0{,}811$ | $-0{,}040$ / $-0{,}042$ | $0{,}385$ / $0{,}464$ | Sem deslocamento em bloco |
| **5. Cognição, Sinapse e Pálio** | 2 / 24 | $-0{,}119$ / $-0{,}120$ | — | $-0{,}016$ / $-0{,}011$ | — | Maioria dos receptores NMDA/AMPA multi-exônicos filtrada pelo limiar estrito de cobertura; no genoma, `CBLN4`, `RTN1`, `NCALD` e `SOX14` lideram |
| **6. Anti-Câncer, Reparo de DNA e Telômeros** | 24 / 25 | $-0{,}020$ / $-0{,}025$ | $0{,}294$ / $0{,}250$ | $-0{,}008$ / $-0{,}013$ | $0{,}800$ / $0{,}771$ | Direções opostas cancelam a média: forte restrição em `PTEN` ($\text{Rho} = -0{,}37\text{ a }-0{,}47, \text{permP} \le 0{,}010$), `RAD50` ($\text{Rho}_{\text{Brain}} = -0{,}20, \text{permP} = 0{,}050$) e `TERF2IP` ($\text{Rho}_{\text{LE}} = -0{,}27, \text{permP} = 0{,}040$ na Smith) vs. aceleração em `NUDT1` ($\text{Rho} = +0{,}27, \text{permP} = 0{,}038$) |

---

## 5. H8 Intra-Clado: Taxa Evolutiva × Longevidade Dentro de Psittaculidae vs. Restante da Ordem

A análise intra-clado (`h8_within_clade_summary.tsv`, `rer_*_h8_{Psittaculidae,rest}_panel.tsv`) complementa o teste entre clados da Fase 2 (HyPhy RELAX e aBSREL):

1. **Dentro de Psittaculidae ($n = 15$ espécies, até 27 ramos, 54 genes do painel com cobertura suficiente):**
   - Com 15 espécies, o poder estatístico intra-clado é limitado e **nenhum gene do painel passa $q < 0{,}10$** (`min panelQ = 0,84`).
   - Em nível nominal ($P \le 0{,}05$ consistente nas duas árvores), destacam-se:
     - **`CAT` (Catalase, M6):** forte restrição purificadora nas espécies mais longevas de Psittaculidae ($\text{Rho} = -0{,}651 / -0{,}640$, $P = 0{,}0085 / 0{,}010$, $\text{permP} = 0{,}058 / 0{,}062$).
     - **`BBC3` (PUMA, M3):** aceleração/relaxamento ($\text{Rho} = +0{,}527 / +0{,}517$, $P = 0{,}044 / 0{,}048$).
     - **`CDKN1B` ($p27^{\text{Kip1}}$, M2):** aceleração ($\text{Rho} = +0{,}393 / +0{,}396$, $P = 0{,}052 / 0{,}050$, $\text{permP} = 0{,}054 / 0{,}044$).
2. **Dentro do restante de Psittaciformes (`rest`, $n = 43$ espécies, até 83 ramos, 58 genes do painel):**
   - No módulo **M8 (Proteostase/Autofagia)** na árvore de Burgio 2019, **2 genes do painel passam `panelQ < 0,10`** (e permanecem $\text{permP} \le 0{,}016$ na árvore de Smith):
     - **`PSMA3` (M8):** restrição purificadora em espécies mais longevas fora de Psittaculidae ($\text{Rho} = -0{,}37 / -0{,}38$, $P \le 0{,}001$, $\text{permP} = 0{,}002 / 0{,}016$, $\text{panelQ} = 0{,}058$ na Burgio).
     - **`PSMA4` (M8):** $\text{Rho} = +0{,}45 / +0{,}47$ dentro de `rest` ($P < 0{,}001$, $\text{permP} = 0{,}002 / 0{,}006$, $\text{panelQ} = 0{,}058$ na Burgio) — notavelmente invertido em relação ao sinal na ordem inteira ($\text{Rho} = -0{,}40$), indicando que o sinal global de `PSMA4` reflete o contraste basal entre Psittaculidae e os demais clados.
   - Além deles, **`PTEN` (M2)** mantém restrição purificadora significativa dentro de `rest` ($\text{Rho} = -0{,}29 / -0{,}32$, $P = 0{,}009 / 0{,}003$, $\text{permP} = 0{,}034$ na Burgio), **`MCL1` (M3)** é conservado ($\text{Rho} = -0{,}32 / -0{,}38$, $P = 0{,}031 / 0{,}010$, $\text{permP} = 0{,}012 / 0{,}008$) e **`MDM4` (M1)** é acelerado ($\text{Rho} = +0{,}23 / +0{,}26$, $P = 0{,}043 / 0{,}018$, $\text{permP} = 0{,}064 / 0{,}018$).

---

## 6. Síntese Integrada do Projeto (Fases 1, 2 e 3) e Conclusões para o Manuscrito

1. **Fase 1 (Auditoria Completa de 126 Montagens e Alometria Filogenética):**
   - Todas as 126 montagens de Psittaciformes (96 Iridian + 30 não-Iridian) foram auditadas com `compleasm` (`aves_odb10`), resultando em **71 espécies incluídas** (`BUSCO >= 70%`, sem reprovação de identidade; 20 Psittaculidae e 51 outros Psittaciformes).
   - Não há viés técnico entre qualidade de montagem (`busco_C`, `contig_N50`, GC, desvio de tamanho) e longevidade relativa ($p \ge 0{,}11$).
   - A redução de longevidade em Psittaculidae é **sensível à topologia filogenética** ($\beta = -0{,}104, p = 0{,}003$ na supertree de Burgio 2019 vs. $\beta = -0{,}076, p = 0{,}050$ na árvore molecular de Smith 2024 com $\lambda = 0{,}88$), pois os Psittaculidae de vida curta concentram-se numa única radiação (Loriini + *Melopsittacus* + *Bolborhynchus*). Portanto, não deve ser apresentada como fato estabelecido independente da árvore.
2. **Fase 2 (Painel de 101 Genes, Nulo Genômico de 298 Genes e Verificação em Reads):**
   - **Sem perda gênica confirmada nem colapso modular (H8 Nível 1 rejeitado):** As 329 chamadas iniciais de genes interrompidos deveram-se quase inteiramente a artefatos de referência única ou fragmentação de contigs curtos/ricos em GC; a checagem com segunda referência (*Strigops habroptila*) e reads Illumina do SRA eliminou todas as perdas fixadas aparentes. Contra o modelo nulo genômico de 298 genes, nenhum módulo M1–M8 sobrevive à correção para múltiplos testes.
   - **Sinal molecular consistente no complexo telomérico (M7) e reparo de quebras duplas (`MRE11`):** Seleção episódica forte e replicável no ramo de *Melopsittacus undulatus* em `RTEL1` ($q = 1{,}5\times 10^{-8}$), `TERF2` ($q = 5\times 10^{-5}$), `TERT` ($q = 0{,}004$), `ACD` ($q = 0{,}005$) e `TERF1` ($q = 0{,}006$), acompanhada de intensificação de seleção em `TERT` ($k = 5{,}5 / 3{,}1$) e `MRE11` ($k = 2{,}6$) em Psittaculidae.
3. **Fase 3 (RERconverge em 7.756 Genes nas Duas Árvores):**
   - **Coevolução genômica entre longevidade e encefalização ($r = 0{,}666\text{--}0{,}667$):** A pressão seletiva associada à longevidade relativa e ao tamanho cerebral relativo atua de forma fortemente acoplada no proteoma dos psitacídeos, mesmo após remover o efeito linear da longevidade sobre o cérebro (`rel_brain_given_le`).
   - **Pleiotropia entre supressão tumoral, proteostase e neurobiologia:** O supressor tumoral mestre **`PTEN`** ($\text{Rho} = -0{,}37\text{ a }-0{,}47$) e o eixo **`CTNNB1`/`CTNNBIP1`** sofrem purificação intensificada nas linhagens mais longevas e encefalizadas (incluindo dentro do clado não-Psittaculidae), enquanto **`MAD2L2` (`REV7`)**, a chaperona de autofagia neuronal **`HSPA8` (`HSC70`)** e organizadores sinápticos/paliais (**`RTN1`**, **`SOX14`**, **`NCALD`**, **`CBLN4`**) exibem aceleração evolutiva coordenada nas duas árvores.
