# Auditoria Red-Team (9 Pontos): Tentando Quebrar os Resultados do Projeto *Parrots*

**Data:** 2026-11-10  
**Repositório / Branch:** `redkk123/Parrots` (`claude/new-session-lwh7xh`)  
**Objetivo:** Submeter as duas principais descobertas preliminares — (A) o acoplamento genômico $r = 0{,}666$ entre longevidade relativa e tamanho cerebral relativo em 7.756 genes (`RERconverge`), incluindo os hits de topo (`PTEN`, `MAD2L2`, `HSPA8`, etc.) e `rel_brain_given_le` ($r = 0{,}639$), e (B) a seleção positiva episódica em genes teloméricos (`M7`) em *Melopsittacus undulatus* (`aBSREL`) — a 9 testes empíricos projetados para invalidar ou falsificar os achados.

---

## Resumo Executivo: O que Sobreviveu vs. O que Quebrou

| Alvo da Auditoria | Resultado Original | Veredito Red-Team | Diagnóstico Empírico Real |
| :--- | :--- | :--- | :--- |
| **1. Reprodutibilidade de $r = 0{,}666$** | $r = 0{,}6663$ | **PASSA (Determinístico)** | `Rho` é 100% determinístico ($\max|\Delta\text{Rho}| = 5{,}0\times 10^{-16}$ com semente `999999` vs `20261028`). Apenas `permP` varia levemente ($\text{cor} = 0{,}991\text{–}0{,}993$). |
| **2. Piso de ruído e controles negativos** | $r = 0{,}6663$ tratado como sinal específico | **QUEBRA PARCIALMENTE** | Traços supostamente não relacionados com resíduo PGLS não branqueado dão $r = +0{,}679$ (massa) e $-0{,}628$ (ninhada bruta), e permulações aleatórias (BM) atingem até $+0{,}670$ / $+0{,}716$ ($p = 0{,}010\text{–}0{,}025$), porque o `RERconverge` amplifica qualquer correlação de ramos ($r_{\text{path}}$) via não-isotropia da matriz de taxas ($\text{cor}(r_{\text{Rho}}, r_{\text{path}}) = 0{,}88\text{–}0{,}91$). Com resíduo OLS de ninhada (`rel_clutch_size`), cai ao piso ($r = -0{,}102$). |
| **3. Viés técnico dos genes** | $r = 0{,}6663$ em 7.756 genes | **PASSA PARA $r$, MAS EXPOSE ARTEFATO NOS HITS** | Regressar 7 a 14 covariáveis técnicas reduz $r$ de $0{,}666$ apenas para $0{,}653$; o sinal persiste em todos os quartis ($0{,}630\text{–}0{,}696$). **Porém**, `total_tree_length` correlaciona com `Rho_Brain` ($r = -0{,}178$, $p = 1{,}8\times 10^{-56}$) e `Rho_LE` ($r = -0{,}154$, $p = 3{,}1\times 10^{-42}$). **Achado crítico:** **89,3% (117/131)** dos genes compartilhados (`P<0,05 & permP<0,05`) têm **>80% dos ramos zerados (`1e-08`)** — incluindo `PTEN` (`99,1%`), `MAD2L2` (`99,0%`), `CTNNB1` (`99,1%`), `HSPA8` (`98,2%`). Neles, o resíduo `RER` espelha principalmente $\pm\sqrt{\text{nv}}$ do genoma de fundo. |
| **4. Dependência de clados** | $r = 0{,}666$ (7.756 genes); $r_{\text{tip}} = 0{,}478$ (50 spp) | **DEPENDE FORTEMENTE DE PSITTACULIDAE** | Tirar qualquer clado isolado (`Strigopidae`, `Cacatuidae`, `Loriini`, `M. undulatus`, `Androglossini`, `Arini`) mantém $r_{\text{Rho}}$ entre $0{,}649$ e $0{,}798$. Mas remover **todos os 13 Psittaculidae** faz a correlação fenotípica cair de $0{,}478$ para **$0{,}195$** (Spearman **$0{,}097$**) e $r_{\text{Rho}}$ cair de $0{,}748$ para **$0{,}522$**. |
| **5. O cérebro residual (`rel_brain_given_le`)** | $r$ cai só de $0{,}666$ para $0{,}639$ | **QUEBRA (METODOLÓGICO)** | Em `scripts/35_prepare_rer_inputs.py` (l. 38–40), usou-se resíduo PGLS ($\lambda = 0{,}90$, $\hat{\beta}_{\text{GLS}} = 0{,}040$) sem branqueamento $V_\lambda^{-1/2}$, deixando `rel_brain_given_le` correlacionado a LE em $r_{\text{tip}} = 0{,}427$ ($\text{cor}(\text{Rho}_{\text{Brain}}, \text{Rho}_{\text{PGLS}}) = 0{,}9987$). Com **ortogonalização OLS/parcial real**, o acoplamento cai pela metade: **$r = 0{,}332$** (OLS nas pontas), **$r = 0{,}389$** (ortogonalização nos ramos) e **$r = 0{,}414$** (correlação parcial analítica). |
| **6. Taxa-base `aBSREL` e correção múltipla** | 5/8 em `M7` e 6/86 no painel | **CONTRASTE COM CONTROLES É REAL ANTES DA LIMPEZA** | BH-FDR em `scripts/32_h8_summary.py` foi corrigido corretamente entre os 86–87 genes do painel (além de Holm-Bonferroni nos 2 ramos). Nos controles, apenas **1/19 (5,3%, $p<0{,}05$)** e **0/19 (0%, $q<0{,}10$)** passaram (`0/298` no RELAX com $\omega>1$). O painel fora de `M7` teve apenas **1/78 (1,3%)** (`MYC`). Só `M7` (telômeros) concentrou **5/8 (62,5%)** ($p_{\text{Fisher}} = 0{,}00069$ em $q<0{,}10$). |
| **7. Alinhamentos e `MultiHit` nos telômeros** | `RTEL1`, `TERF2`, `TERT`, `ACD`, `TERF1` significativos em *M. undulatus* | **QUEBRA PARA 3 GENES; TODOS ZERAM COM MULTI-HIT** | Realinhar com **MAFFT L-INS-i** mantém o sinal. Porém: (a) **limpar janelas ruidosas/gaps** elimina `RTEL1` ($1{,}7\times 10^{-10} \to 1{,}00$), `ACD` ($2{,}5\times 10^{-4} \to 0{,}125$) e `TERF2` ($1{,}8\times 10^{-6} \to 0{,}068$), restando apenas `TERT` ($0{,}0052$) e `TERF1` ($0{,}0037$); (b) ativar **`--multiple-hits Double+Triple`** zera **todos** (`RTEL1`=0,967; `TERF2`=1,00; `TERT`=1,00; `TERF1`=1,00; `ACD`=0,050/0,132). Ortologia dupla (periquito × kakapo) confirma ortologia real em 100% (`RTEL1` 40/40, `TERF2` 30/30, `TERT` 35/35, `TERF1` 27/27, `POT1` 27/27; `ACD` 8/11 em 2 loci). |
| **8. Replicação nas 13 novas espécies (Holdout)** | Painel de 101 genes em 71 spp + 30 splits CV | **REPLICAM NO PAINEL (`r=0,94`), NÃO ISOLADOS NAS 13 SPP** | As 13 espécies novas só têm extrações do painel de 97 genes (não dos 7.756 genes genômicos). No painel em 71 spp, `Rho_LE` replica muito bem (`cor = 0,9416` vs 58 spp; `TERT`, `TERF2`, `RTEL1`, `TERF1`, `POT1` todos negativos e significativos em 71 spp), mas isoladas nas 13 espécies novas ($n = 13$, ~10–24 ramos) quase nenhum gene atinge $p<0{,}05$ por baixo poder. **`PTEN` inverte de sinal** (`-0,374` na Fase 3 genômica vs `+0,184` no painel de 91 genes por depender da normalização de fundo `nv` em um gene com 99% de ramos zerados). Nos 30 splits CV dentro de 50 spp, `PTEN`, `MAD2L2`, `HSPA8` mantêm sinal em 97–100% porque compartilham o mesmo vetor `RER` pré-computado. |
| **9. Honestidade estatística** | Inventário completo de todos os testes | **78.802 TESTES DOCUMENTADOS** | Inventário completo em [`9_statistical_honesty_inventory.tsv`](file:///home/red/parrots/results/redteam_audit/9_statistical_honesty_inventory.tsv) incluindo todos os achados nulos (dieta, latitude, qualidade, pseudogenização H8a, stem Psittaculidae, módulos M1–M8 contra o nulo, FDR genômico de `permQ`, e H8 intra-clado). |

---

## Resultados Detalhados por Questão

### 1. Reprodutibilidade (`r = 0,666`)
- **Scripts exatos:** [`scripts/35_prepare_rer_inputs.py`](file:///home/red/parrots/scripts/35_prepare_rer_inputs.py) (prepara fenótipos e árvores mestre) e [`scripts/R/rerconverge_h6_h8.R`](file:///home/red/parrots/scripts/R/rerconverge_h6_h8.R) (`set.seed(20261028)`, l. 63).
- **Arquivos de saída originais:** [`rer_burgio2019_rel_log_life_expectancy.tsv`](file:///home/red/parrots/results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv) e [`rer_burgio2019_rel_log_brain.tsv`](file:///home/red/parrots/results/phase3_genome/rer_burgio2019_rel_log_brain.tsv).
- **Re-execução com semente diferente (`999999`)** via [`scripts/R/redteam_rer_audit.R`](file:///home/red/parrots/scripts/R/redteam_rer_audit.R) → [`1_reproducibility.tsv`](file:///home/red/parrots/results/redteam_audit/1_reproducibility.tsv):
  - `r(Rho_LE, Rho_Brain)` salvo (`20261028`): **`0,666303820424148`**
  - `r(Rho_LE, Rho_Brain)` com semente `999999`: **`0,666303820424147`** ($\max|\Delta\text{Rho}| = 4{,}996\times 10^{-16}$, zero de máquina).
  - Na matriz `RER` idêntica das 50 espécies com ambos os traços (`rers_50`): **`r = 0,7478`**.
  - A semente aleatória só afeta as **permulações** (`permP`): $\text{cor}(\text{permP}_{\text{LE}}) = 0{,}9907$ (diferença mediana $|\Delta| = 0{,}0259$), $\text{cor}(\text{permP}_{\text{Brain}}) = 0{,}9932$ (mediana $|\Delta| = 0{,}0200$). Genes compartilhados (`P<0,05 & permP<0,05`, mesmo sinal): `131` (semente `20261028`) vs `135` (semente `999999`), com `124` em comum (`94,7%`).

---

### 2. Piso de Ruído do Acoplamento e Controles Negativos
Arquivo: [`2_coupling_null_and_controls.tsv`](file:///home/red/parrots/results/redteam_audit/2_coupling_null_and_controls.tsv).

| Comparação / Modelo Nulo | $r_{\text{pontas}}$ | $r_{\text{ramos}}$ ($r_{\text{path}}$) | $r(\text{Rho}_1, \text{Rho}_2)$ genômico | Nota |
| :--- | :---: | :---: | :---: | :--- |
| **Observado:** `rel_log_life_expectancy` × `rel_log_brain` (salvo) | +0,4783 | +0,4154 | **+0,6663** | Matriz RER própria por traço |
| **Observado (50 spp pareado):** `LE` × `Brain` | +0,4783 | +0,4253 | **+0,7478** | Mesma matriz `rers_50` |
| **Controle 1:** `rel_log_life_expectancy` × `log_body_mass` ($n=58$) | +0,4102 | +0,2885 | **+0,6792** | Resíduo PGLS não branqueado retém efeito de massa! |
| **Controle 2:** `rel_log_life_expectancy` × `clutch_size` bruta ($n=51$) | -0,4345 | -0,3482 | **-0,6280** | Ninhada bruta é quase espelho negativo |
| **Controle 3:** `rel_log_life_expectancy` × `rel_clutch_size` (resíduo OLS na massa, $n=51$) | -0,1725 | -0,2145 | **-0,1024** | Quando ortogonalizado à massa, vai ao piso! |
| **Controle 4:** `rel_log_life_expectancy` × `rel_incubation_days` (resíduo OLS na massa, $n=50$) | +0,2060 | +0,0513 | **+0,3417** | Correlação residual moderada |
| **Nulo 2A:** 200 pares de permulações BM independentes | ~0 | média -0,002 | média -0,009 (DP=0,2885), **IC 95% [-0,535; +0,529]**, **máx +0,6703** | **$p_{\text{emp}}(r \ge 0{,}666) = 0{,}010$ (2/201)**; $\text{cor}(r_{\text{Rho}}, r_{\text{path}}) = 0{,}8831$ |
| **Nulo 2B:** 200 pares simulados juntos (BM acoplado empírico) | -0,2786 | -0,2762 | **IC 95% [-0,791; +0,453]**, **máx +0,6601** | $\text{cor}(r_{\text{Rho}}, r_{\text{path}}) = 0{,}9133$ |
| **Nulo 2C:** `LE` observado fixo × 200 permulações de `Brain` | ~0 | média +0,033 | média +0,019 (DP=0,3504), **IC 95% [-0,594; +0,651]**, **máx +0,7155** | **$p_{\text{emp}}(r \ge 0{,}666) = 0{,}0249$ (5/201)** |

**Diagnóstico:** O valor $0{,}666$ fica na cauda superior da distribuição nula ($p_{\text{emp}} = 0{,}010$ a $0{,}025$), mas a dispersão sob o nulo é enorme ($\text{DP} = 0{,}29\text{–}0{,}35$, podendo chegar a $0{,}716$ por puro acaso). O motivo é estrutural: quaisquer dois traços que compartilhem correlação nos ramos da árvore ($r_{\text{path}}$) têm essa correlação amplificada para $r_{\text{Rho}}$ ($\text{cor} = 0{,}88\text{–}0{,}91$).

---

### 3. Viés Técnico dos Genes (e o Artefato dos Genes Invariantes)
Arquivos: [`3_technical_confounders.tsv`](file:///home/red/parrots/results/redteam_audit/3_technical_confounders.tsv) e [`3_coupling_after_confounder_removal.tsv`](file:///home/red/parrots/results/redteam_audit/3_coupling_after_confounder_removal.tsv) (gerados por [`45_redteam_python_audit.py`](file:///home/red/parrots/scripts/45_redteam_python_audit.py)).

- **Correlações técnicas nos 7.756 genes:**
  - Comprimento (`aln_codons`): $r = -0{,}0921$ com $|\text{Rho}_{\text{LE}}|$ ($p = 4{,}5\times 10^{-16}$), mas $r = +0{,}0023$ com $\text{Rho}_{\text{LE}}$ com sinal ($p = 0{,}84$).
  - Taxa evolutiva de fundo (`total_tree_length`): **$r = -0{,}1537$** com $\text{Rho}_{\text{LE}}$ ($p = 3{,}1\times 10^{-42}$) e **$r = -0{,}1784$** com $\text{Rho}_{\text{Brain}}$ ($p = 1{,}8\times 10^{-56}$); $r = -0{,}1179$ com $|\text{Rho}_{\text{LE}}|$.
  - Conteúdo GC / GC3: $r = -0{,}0479$ com $\text{Rho}_{\text{LE}}$ ($p = 2{,}5\times 10^{-5}$).
  - Fração de espécies faltando (`missing_fraction`): $r = +0{,}0658$ com $|\text{Rho}_{\text{LE}}|$ ($p = 6{,}7\times 10^{-9}$) e $r = +0{,}1005$ com $|\text{Rho}_{\text{Brain}}|$ ($p = 7{,}5\times 10^{-19}$).
  - Fração de gaps: $r = -0{,}0575$ com $\text{Rho}_{\text{LE}}$ ($p = 4{,}0\times 10^{-7}$).
- **Acoplamento após remover efeitos técnicos:**
  - Bruto: **$0{,}6663$** (Spearman $0{,}6323$; $r(|\text{Rho}|) = 0{,}5430$)
  - Após regressão OLS de 7 covariáveis lineares + log: **$0{,}6539$**
  - Após regressão OLS quadrática (14 covariáveis): **$0{,}6532$**
  - Dentro de todos os quartis de comprimento, taxa de fundo, GC3 e espécies faltantes: **$r = 0{,}6297\text{ a }0{,}6965$**.
- **DESCOBERTA CRÍTICA DA AUDITORIA (Genes com Ramos Zerados):**
  - Quando estratificamos os 7.756 genes pela fração de ramos com comprimento nulo (`<= 1e-6`, isto é, `1e-08` no LG):
    - **Todos os 7.756 genes:** $r(\text{Rho}) = 0{,}6663$, **131** hits compartilhados (`P < 0,05 & permP < 0,05`).
    - **Ultraconservados (>80% ramos zerados, $n=1.450$):** $r = 0{,}6761$, **117/131 hits compartilhados (`89,3%`!)**. Exemplos do Top 20: `PTEN` (**110/111 ramos zerados, 99,1%**), `MAD2L2` (**104/105, 99,0%**), `HSPA8` (**107/109, 98,2%**), `CTNNB1` (**112/113, 99,1%**), `PSMA4` (**110/111, 99,1%**), `SOX14` (**98,1%**), `RTN1` (**97,3%**).
    - **Genes informativos ($\ge 70\%$ ramos não-nulos, $n=1.263$):** $r = 0{,}6982$, mas apenas **2 genes** (`0,16%`) passam `P < 0,05 & permP < 0,05` nos dois traços! (Com $\ge 20$ ramos não-nulos, $n=6.182$: apenas **14 genes**, encabeçados por `TKT`, `PHOSPHO1`, `TADA2A`, `SLC4A7`, `HUS1`).
  - *Implicação:* O acoplamento global ($0{,}65\text{–}0{,}69$) não é criado por viés de comprimento/GC, mas **quase 90% dos genes de topo compartilhados (`PTEN`, `MAD2L2`, `HSPA8`, etc.) são genes invariantes com ~98–99% dos ramos zerados**, cujo resíduo `RER` deriva matematicamente da subtração do vetor médio genômico $\sqrt{\text{nv}}$!

---

### 4. Quem Puxa o Sinal (Análise Leave-One-Clade-Out)
Arquivo: [`4_leave_one_clade_out.tsv`](file:///home/red/parrots/results/redteam_audit/4_leave_one_clade_out.tsv).

- **Correlação entre `rel_log_life_expectancy` e `rel_log_brain` nas 50 espécies:** Pearson **$r = 0{,}4783$** ($r_{\text{Spearman}} = 0{,}4767$, $r_{\text{path}} = 0{,}4253$, $r_{\text{Rho},50\text{spp}} = 0{,}7478$).
- **Removendo um clado por vez:**
  - Tirar `Strigopidae` (-2 spp): $r_{\text{tip}} = 0{,}5095$ | $r_{\text{Rho}} = 0{,}7646$
  - Tirar `Cacatuidae` (-6 spp): $r_{\text{tip}} = 0{,}4631$ | $r_{\text{Rho}} = 0{,}7626$
  - Tirar `Loriini` (-5 spp): $r_{\text{tip}} = 0{,}4773$ | $r_{\text{Rho}} = 0{,}7483$
  - Tirar `Melopsittacus undulatus` (-1 sp): $r_{\text{tip}} = 0{,}4458$ | $r_{\text{Rho}} = 0{,}7285$
  - Tirar `Androglossini` (-9 spp): $r_{\text{tip}} = 0{,}5419$ | $r_{\text{Rho}} = 0{,}7975$
  - Tirar `Arini` (-12 spp): $r_{\text{tip}} = 0{,}3798$ | $r_{\text{path}} = 0{,}2407$ | $r_{\text{Rho}} = 0{,}6489$
  - **Tirar TODOS os `Psittaculidae` (-13 spp, restando 37 espécies):** $r_{\text{tip}}$ desaba de **$0{,}4783$ para $0{,}1955$** (Spearman **$0{,}0973$**!), $r_{\text{path}}$ cai para **$0{,}2095$**, e $r_{\text{Rho}}$ cai de **$0{,}7478$ para $0{,}5225$**. (Dentro das 37 não-Psittaculidae, `PTEN` e `MAD2L2` ainda mantêm `Rho = -0,356` e `+0,344`).
- *Conclusão:* Nenhum gênero sozinho (como *Melopsittacus*) sustenta o acoplamento, mas **o contraste entre Psittaculidae (vida mais curta / cérebro relativamente menor) e o restante da ordem** responde pela maior parte da associação fenotípica linear entre cérebro e longevidade.

---

### 5. O Cérebro Residual (`rel_brain_given_le`) e Correlação Parcial Real
Arquivos: [`5_residual_brain_summary.tsv`](file:///home/red/parrots/results/redteam_audit/5_residual_brain_summary.tsv) e [`5_residual_brain_key_genes.tsv`](file:///home/red/parrots/results/redteam_audit/5_residual_brain_key_genes.tsv).

- **Como foi calculado originalmente ([`scripts/35_prepare_rer_inputs.py`](file:///home/red/parrots/scripts/35_prepare_rer_inputs.py), l. 38–40):**
  `p0b.subset_fit(full, both, "rel_log_brain", ["rel_log_life_expectancy"])` ajustou um **PGLS com Pagel's $\lambda = 0{,}90$**, cuja inclinação filogenética foi de apenas **$\hat{\beta}_{\text{GLS}} = 0{,}040 \pm 0{,}032$**, e subtraiu `rel_log_brain - 0.040 * rel_log_life_expectancy` **sem transformar/branquear por $V_\lambda^{-1/2}$**.
- **Por que o $r$ só caiu de $0{,}666$ para $0{,}639$?**
  Porque subtrair `0,040 * LE` quase não mudou o vetor fenotípico: `rel_brain_given_le` continuou correlacionado com `rel_log_life_expectancy` em **$r_{\text{tip}} = 0{,}4269$** ($r_{\text{path}} = 0{,}3704$), e $\text{cor}(\text{Rho}_{\text{Brain}}, \text{Rho}_{\text{Brain|LE,PGLS}}) = \mathbf{0{,}9987}$! O teste antigo era essencialmente redundante.
- **Correlação parcial / ortogonalização de verdade:**
  1. **Resíduo OLS nas pontas** (`rel_brain_given_le_OLS`, $\hat{\beta}_{\text{OLS}} = 0{,}258$, $r_{\text{tip}} = 0{,}000$, $r_{\text{path}} = -0{,}063$): $r(\text{Rho}_{\text{LE}}, \text{Rho}_{\text{Brain|LE,OLS}})$ cai pela metade para **$0{,}3322$** ($0{,}3691$ na matriz pareada de 50 spp).
  2. **Ortogonalização no nível dos ramos** (`paths_brain_orth_le`, $r_{\text{path}} = 0{,}000$ exato): cai para **$0{,}3888$** ($0{,}4332$ em 50 spp).
  3. **Correlação parcial analítica por gene** $r(\text{RER}_g, P_{\text{Brain}} \mid P_{\text{LE}})$: **$0{,}4135$** ($0{,}4612$ em 50 spp).
  4. *Nos genes-chave (`5_residual_brain_key_genes.tsv`):* `MAD2L2` ainda retém sinal paramétrico e de permulação para `Brain_OLS` ($\text{Rho} = +0{,}455$, $P = 5{,}75\times 10^{-6}$, $\text{permP} = 0{,}0297$); `ANAPC11` tem $\text{permP} = 0{,}0396$; já `PTEN` ($\text{permP} = 0{,}109$), `CTNNB1` ($\text{permP} = 0{,}228$) e `HSPA8` ($\text{permP} = 0{,}198$) perdem significância na permulação após a ortogonalização real.

---

### 6. Taxa-Base da Seleção Positiva (`aBSREL`) nos Controles e Correção Múltipla
Arquivos: [`6_absrel_base_rate.tsv`](file:///home/red/parrots/results/redteam_audit/6_absrel_base_rate.tsv) e [`6_absrel_null_genes_detail.tsv`](file:///home/red/parrots/results/redteam_audit/6_absrel_null_genes_detail.tsv).

- **Como o $q$-value foi corrigido:** Verificado em [`scripts/32_h8_summary.py`](file:///home/red/parrots/scripts/32_h8_summary.py) (l. 38–46 e 75–77). Cada arquivo `.absrel.json` do HyPhy já traz o $P$ com **correção Holm-Bonferroni entre os 2 ramos testados** por gene (`Melopsittacus_undulatus` e o ramo-tronco de `Psittaculidae`), e o script `32_h8_summary.py` aplicou **BH-FDR entre todos os 86–87 genes do painel** (`groupby("tree")[col].transform(bh)`), não apenas entre os 2 ramos.
- **Taxa-base nos genes de controle (`hyphy_null/burgio2019`):**
  - Nos **19 genes de controle** executados no `aBSREL`: apenas **1/19 (5,3%)** teve `p_holm2 < 0,05` em *M. undulatus* (`NDUFC1`, $p = 0{,}0385$) e **0/19 (0,0%)** sobreviveu ao BH-FDR ($q = 0{,}626$).
  - Além disso, em **todos os 298 genes de controle** já rodados no RELAX, **0/298 (0,0%)** apresentaram $\omega > 1$ nas pontas de `Psittaculidae`.
  - **Comparação com o painel H8 (86 genes testáveis):** No painel completo, **6/86 (7,0%)** passaram ($p_{\text{Fisher}} = 1{,}0$ vs controles). Fora dos telômeros (`M1–M6, M8`), apenas **1/78 (1,3%)** passou (`MYC`). **Somente o módulo telomérico `M7` (8 genes)** teve **5/8 (62,5%)** com $q < 0{,}10$ ($p_{\text{Fisher}} = 0{,}00369$ em $p < 0{,}05$; $p_{\text{Fisher}} = 0{,}000694$ em $q < 0{,}10$).

---

### 7. Auditoria dos Alinhamentos Teloméricos (`RTEL1`, `TERF2`, `TERT`, `ACD`, `TERF1`) e Ortologia
Arquivos: [`7_telomere_alignment_and_rbh.tsv`](file:///home/red/parrots/results/redteam_audit/7_telomere_alignment_and_rbh.tsv) e [`7_telomere_substitution_windows.txt`](file:///home/red/parrots/results/redteam_audit/7_telomere_substitution_windows.txt).

| Gene | $n_{\text{sp}}$ | Subst. AA *M. undulatus* (Isoladas / Janela Ruidosa ou Gap) | Códons Mascarados | `aBSREL p` Original (MAFFT) | `aBSREL p` (MAFFT L-INS-i) | `aBSREL p` (Limpeza de Janelas) | `aBSREL p` (`--multiple-hits Double+Triple`) | `aBSREL p` (Limpo + Multi-Hit) | Ident. *Melo* RefSeq / N loci | Concordância Dupla (*Kakapo*) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **RTEL1** | 23 | 279 (58 / **221**) | 587 / 1307 | $1{,}73\times 10^{-10}$ | $1{,}28\times 10^{-9}$ | **1,0000** | **0,9670** | **1,0000** | 1,000 / 2 loci | **40/40** |
| **TERF2** | 18 | 89 (32 / **57**) | 299 / 646 | $1{,}79\times 10^{-6}$ | $2{,}44\times 10^{-6}$ | **0,0677** | **1,0000** | **1,0000** | 1,000 / 1 locus | **30/30** |
| **TERT** | 15 | 196 (95 / 101) | 516 / 1325 | $2{,}63\times 10^{-4}$ | $1{,}12\times 10^{-4}$ | **0,00516** | **1,0000** | **1,0000** | 0,990 / 1 locus | **35/35** |
| **ACD** | 6 | 102 (42 / **60**) | 178 / 658 | $2{,}46\times 10^{-4}$ | $1{,}23\times 10^{-4}$ | **0,1248** | **0,0502** | **0,1320** | 1,000 / 2 loci | **8/11** |
| **TERF1** | 9 | 23 (15 / 8) | 19 / 333 | $3{,}98\times 10^{-4}$ | $1{,}99\times 10^{-4}$ | **0,00366** | **1,0000** | **0,3151** | 1,000 / 1 locus | **27/27** |
| **POT1** | 8 | 84 (52 / 32) | 144 / 780 | 1,0000 | 1,0000 | 1,0000 | 1,0000 | 1,0000 | 1,000 / 1 locus | **27/27** |

**Conclusões desta auditoria:**
1. **Trocar apenas o alinhador (MAFFT L-INS-i) não derruba o sinal:** todos os 5 genes permanecem significativos com L-INS-i.
2. **Inspeção visual (`7_telomere_substitution_windows.txt`) e limpeza de janelas ruidosas:**
   - Em `TERF2` (códons 9–54) e `ACD` (códons 1–63, rico em prolina/glicina no éxon 1 com gaps enormes), as substituições caem em **regiões hipervariáveis e cheias de indels**. Ao limpar colunas com $>15\%$ de gaps e janelas de 10 códons com $\ge 3$ diferenças de aminoácido, **`RTEL1` ($p = 1{,}00$), `ACD` ($p = 0{,}125$) e `TERF2` ($p = 0{,}068$) perdem significância**.
   - Apenas **`TERT` ($p = 0{,}0052$)** e **`TERF1` ($p = 0{,}0037$)** resistem à limpeza de janelas!
3. **Modelo com múltiplas substituições por códon (`--multiple-hits Double+Triple`):**
   - **Todos os 5 genes desabam para $p = 1{,}00$ (ou `0,0502` / `0,132` no `ACD`)!** Isso demonstra que a verossimilhança do `aBSREL` padrão foi puxada por códons com **2 ou 3 substituições nucleotídicas simultâneas** (ex. `TERF2` códon 15 `AGG→TCT`, `TERT` códon 14 `ACG→GCA`, `TERF1` códon 247 `CGT/TGT→CAT`), que o modelo padrão de códon único interpreta como $\omega \gg 1$ episódico, mas que o modelo Multi-Hit acomoda sem seleção positiva.
4. **Ortologia (RBH / Referência dupla com *Strigops habroptila*):** A identidade e o mapeamento por coordenadas genômicas (`kakapo` vs `periquito`) confirmam **100% de concordância de lócus** para `RTEL1` (40/40), `TERF2` (30/30), `TERT` (35/35), `TERF1` (27/27) e `POT1` (27/27) — ou seja, não é viés de parálogo causado por usar o proteoma do periquito (exceto `ACD`, com apenas 8/11 concordantes e 2 loci no periquito).

---

### 8. Replicação nas 13 Novas Espécies (Holdout) e Validação Cruzada
Arquivos: [`8_holdout_13_species_list.tsv`](file:///home/red/parrots/results/redteam_audit/8_holdout_13_species_list.tsv), [`8_holdout13_and_71sp_panel_replication.tsv`](file:///home/red/parrots/results/redteam_audit/8_holdout13_and_71sp_panel_replication.tsv) e [`8_holdout_cv_key_genes.tsv`](file:///home/red/parrots/results/redteam_audit/8_holdout_cv_key_genes.tsv) (gerados por [`scripts/R/redteam_replication_71sp.R`](file:///home/red/parrots/scripts/R/redteam_replication_71sp.R)).

- **Plano congelado das 13 novas espécies de teste (`holdout_13`):** 5 `Psittaculidae` (*Agapornis nigrigenis*, *A. fischeri*, *A. roseicollis*, *Neophema chrysogaster*, *Psittacula krameri*) e 8 `rest` (*Psittrichas fulgidus*, *Amazona guildingii*, *A. vittata*, *A. ochrocephala*, *Psittacara finschi*, *P. erythrogenys*, *Forpus passerinus*, *Probosciger aterrimus*).
- **O que existe extraído nessas 13 espécies:** Apenas os **97 genes do painel H8** (`results/phase2_genes_71sp/aa/`), **não** os 7.756 genes do genoma completo (portanto `MAD2L2`, `HSPA8`, `CTNNB1`, `RTN1`, etc. ainda não têm árvores nas 13 novas espécies; `PTEN` e `PSMA4` estão no painel).
- **Resultado no painel de 91 genes válidos (58 spp vs Holdout 13 vs 71 spp):**
  - Concordância global do painel de 58 spp para 71 spp: **$\text{cor}(\text{Rho}_{\text{LE},58}, \text{Rho}_{\text{LE},71}) = 0{,}9416$** e **$\text{cor}(\text{Rho}_{\text{Brain},58}, \text{Rho}_{\text{Brain},71}) = 0{,}9131$**.
  - Nas 13 espécies isoladas, o $N$ de ramos é pequeno (10–24 ramos): apenas 2/91 genes atingem $p < 0{,}05$ para LE (`CDKN1B` $\text{Rho} = -0{,}591$, $p = 0{,}0077$; `BECN1` $+0{,}507$, $p = 0{,}027$; `BNIP3` $-0{,}444$, $p = 0{,}030$).
  - Teloméricos/reparo no painel em 71 spp mantêm direção negativa forte para LE: `RTEL1` (`Rho = -0,589`, $p = 1{,}6\times 10^{-7}$; holdout 13 = `-0,431`, $p = 0{,}095$), `TERT` (`-0,509`, $p = 1{,}5\times 10^{-5}$; holdout 13 = `-0,341`), `TERF2` (`-0,501`, $p = 9{,}7\times 10^{-5}$; holdout 13 = `-0,206`), `TERF1` (`-0,468`, $p = 9{,}1\times 10^{-4}$; holdout 13 = `-0,255`), `POT1` (`-0,430`, $p = 0{,}0020$; holdout 13 = `-0,374`, $p = 0{,}095$).
  - **`PTEN` e `PSMA4` no teste real:** `PSMA4` manteve sinal negativo em 71 spp (`Rho_LE_58 = -0,277`, `Rho_LE_71 = -0,209`, $p = 0{,}014$), embora nas 13 holdout isoladas tenha ficado em `+0,030` ($p = 0{,}889$). Já **`PTEN` mudou de sinal** (`-0,374` no genoma completo de 7.756 genes vs **`+0,184`** no painel de 91 genes em 58 spp, **`-0,206`** $p=0{,}357$ nas 13 holdout e **`+0,140`** $p=0{,}103$ nas 71 spp) — o que expõe a fragilidade apontada na Questão 3 (como `PTEN` tem **110/111 ramos zerados**, seu `RER` depende inteiramente do vetor de normalização `nv` de fundo!).
- **Validação cruzada estratificada (30 splits de 35 treino / 15 teste nas 50 espécies da Fase 3, 7.756 genes):**
  Quando mantemos a mesma matriz `RER` das 50 espécies (`8_holdout_cv_key_genes.tsv`), `PTEN` (`-0,369`, concordância de sinal `100%`), `MAD2L2` (`+0,471`, `100%`), `HSPA8` (`+0,381`, `100%`), `CTNNB1` (`-0,434`, `100%`), `RTN1` (`+0,459`, `100%`) e `SOX14` (`+0,466`, `100%`) replicam em quase todos os subconjuntos de espécies — mas isso mede estabilidade às espécies **dentro da mesma estimativa de `nv`**, não independência do artefato de ramos zerados.

---

### 9. Honestidade Estatística (Inventário Completo)
Arquivo: [`9_statistical_honesty_inventory.tsv`](file:///home/red/parrots/results/redteam_audit/9_statistical_honesty_inventory.tsv).

- **Total aproximado de testes realizados no projeto:** **78.802 testes** (48 na Fase 0/0b + 34 na Fase 1 + 348 nos pilotos H8 de 12/17 spp + 812 na Fase 2 principal de 58 spp + 77.560 na Fase 3 `RERconverge`).
- **Tudo o que foi testado e deu NULO / negativo:**
  1. Efeito da dieta sobre longevidade residual ($p > 0{,}25$);
  2. Latitude / insularidade sobre longevidade residual;
  3. Regimes OU múltiplos fora da divisão principal;
  4. Qualidade de montagem (`N50`, BUSCO) vs longevidade residual ($r = 0{,}04$, $p = 0{,}76$);
  5. Taxa de substituição mitocondrial ($p > 0{,}30$) e distância k-mer Mash genômica após controle da altura da árvore;
  6. **H8a (perda gênica / pseudogenização em Psittaculidae):** **NULO** (100% dos 329 chamados "disrupted" eram artefatos de montagem, polimorfismo heterozigoto como `PSMD6`, ou ausência de cobertura em GC extremo como `ACD`);
  7. **H8b (relaxamento modular no RELAX contra o nulo de 298 genes):** **NULO** — heterogeneidade de $\log_2(k)$ foi nula ($p = 0{,}164$ em Burgio, $p = 0{,}955$ em Smith); heterogeneidade de `prop_relaxed` ($p = 0{,}023$ Burgio) falhou na Smith ($p = 0{,}127$); `M2` ($p = 0{,}024$) não sobreviveu a correção para 16 testes nem à sensibilidade pela mediana ($p = 0{,}07$); módulos `M1, M3–M8` todos nulos ($p = 0{,}09\text{–}1{,}00$);
  8. **H8c (`aBSREL` no ramo-tronco de Psittaculidae):** **NULO** (**0/87** genes com $q < 0{,}10$);
  9. **Fase 3 `permQ < 0,10` para genes individuais nos 7.756 genes:** **NULO** (com 1.000 permulações, `min(permP) = 1/1001 = 0,000999`, insuficiente para sobreviver ao BH-FDR de 7.756 genes);
  10. **H8 intra-clado (`RERconverge` nos 86 genes do painel dentro de `Psittaculidae` e dentro de `rest`):** **NULO** após BH de permulação (**0/86** genes com `panelQ < 0,10`).
