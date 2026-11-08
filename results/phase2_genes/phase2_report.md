# Análise Confirmatória da H8: 58 Espécies Auditadas (8 nov. 2026)

**Análise principal confirmatória** (plano pré-registrado, Histórico de 28 out.). São 58 espécies auditadas com genoma único por espécie (BUSCO ≥ 70 %, identidade mitocondrial e nuclear validadas, preferência por montagens cromossômicas de referência). Inclui 15 espécies de Psittaculidae (incluindo periquito-australiano *Melopsittacus undulatus*, lóris *Lorius lory*, *L. garrulus*, *L. chlorocercus*, *Trichoglossus*, agapornis *Agapornis lilianae*, *A. personatus*, *A. canus*, *Loriculus galgulus*, *Cyanoramphus*, *Lathamus*, *Platycercus*, *Psephotus*, *Neopsephotus*, *Psittacula derbiana*) contra 43 outras espécies de Psittaciformes e outgroup.

Foram avaliados os 101 genes pré-congelados dos módulos M1–M8 (plano §10) nas duas topologias mestras: **Burgio et al. 2019** (supertree) e **Smith et al. 2024** (árvore filogenômica de sensibilidade).

---

## 1. RELAX: Tronco e Ramos Internos de Psittaculidae vs. Referência (FDR q < 0,10)

Resultados consistentes nas duas árvores filogenéticas:

| Direção da Seleção | Genes no Painel (FDR q < 0,10) | Módulos e Valores Representativos |
|---|---|---|
| **Relaxados** (\(k < 1\)) | **TSC2** (\(k \approx 0\), \(q = 0,018\)), **DKC1** (\(k = 0,49\), \(q = 0,004\)), **RTEL1** (\(k \to 0\), \(q = 0,005\)), **NQO1** (\(k = 0,33\), \(q = 0,049\)), **BBC3** (PUMA; \(k \approx 0\), \(q = 0,018\)) | M2 (supressores), M7 (telômeros), M5 (estresse), M3 (apoptose) |
| **Intensificados** (\(k > 1\)) | **TERT** (\(k = 3,1\) a \(5,5\), \(q = 0,005\)), **MRE11** (\(k = 2,6\), \(q = 0,022\)), **SOD1** (\(k = 31,6\) a \(50,0\), \(q = 0,008\)) | M7 (telômeros), M4 (reparo HR), M6 (antioxidante) |

*(Nota: CCNE1 foi intensificado na árvore de Smith e relaxado na de Burgio; apenas os genes acima foram concordantes em magnitude e significância em ambas).*

---

## 2. Comparação com o Nulo Genômico (298 Genes do Universo)

O modelo nulo avaliou o RELAX em 298 genes controle do universo genômico sob a mesma parametrização (`scripts/37_h8_null.py`):

| Módulo | Genes no Painel | Proporção Relaxada | Média Nulo Relax. | p (Proporção) | \(\text{mean}(\log_2 k)\) Painel | \(\text{mean}(\log_2 k)\) Nulo | p (\(\log_2 k\)) |
|---|---|---|---|---|---|---|---|
| **M1** (Checkpoints DNA) | 9 | 0,00 | 0,07 | 1,000 | -0,10 | -0,71 | 0,593 |
| **M2** (Crescimento/Supressores) | 6 | 0,33 | 0,06 | 0,092 | **-3,69** | **-0,49** | **0,024\*** |
| **M3** (Apoptose) | 10 | 0,10 | 0,04 | 0,663 | -1,37 | -0,00 | 0,094 |
| **M4** (Reparo DNA) | 11 | 0,00 | 0,09 | 0,805 | -0,05 | -0,79 | 0,378 |
| **M5** (Imunidade/Estresse) | 8 | 0,12 | 0,06 | 0,803 | +0,35 | -0,68 | 0,188 |
| **M6** (Metabolismo/Antiox.) | 13 | 0,00 | 0,03 | 1,000 | -0,49 | -0,30 | 0,713 |
| **M7** (Telômeros) | 8 | 0,25 | 0,08 | 0,328 | -1,13 | -0,74 | 0,641 |
| **M8** (Autofagia/Proteostase) | 22 | 0,00 | 0,03 | 1,000 | -0,27 | +0,04 | 0,585 |

- **M2 (Crescimento / Supressores de Tumor)** apresenta relaxamento significativamente mais intenso que a taxa de fundo do genoma (\(p = 0,024\)), puxado especialmente pela via mTOR/supressão (*TSC2*).
- **M3 (Apoptose)** exibe tendência marginal ao relaxamento médio (\(p = 0,094\)), puxado por *BBC3* (PUMA).
- Não há colapso generalizado nem perda indiscriminada em todos os módulos: M1, M4, M6 e M8 mantêm seleção purificadora idêntica ou superior ao nulo genômico.

---

## 3. aBSREL: Seleção Episódica (q < 0,10)

- **Ramo Terminal do Periquito (*Melopsittacus undulatus*):**
  - **RTEL1** (\(q = 1,5 \times 10^{-8}\))
  - **MYC** (\(q = 7,3 \times 10^{-7}\))
  - **TERF2** (\(q = 5,2 \times 10^{-5}\))
  - **TERF1** (\(q = 5,8 \times 10^{-3}\))
  - **TERT** (\(q = 4,6 \times 10^{-3}\))
  - **ACD** (\(q = 4,6 \times 10^{-3}\))
- **Tronco de Psittaculidae:**
  - Nenhum gene isolado atingiu significância estrita com \(q < 0,10\) no tronco isolado das 58 espécies; os sinais concentram-se no clado como um todo (RELAX) ou na aceleração específica do ramo terminal de vida curta (*Melopsittacus*).

---

## 4. Auditoria de Integridade e Checagem com Reads (*Read-Check*)

A triagem inicial identificou 31 chamadas de genes interrompidos (*disrupted*) concordantes entre referências (periquito e kakapo), concentradas em 6 candidatos a pseudogenização com cobertura proteica \(\ge 90\%\) e \(\le 2\) stops. Todas foram inspecionadas com reads brutos Illumina independentes (SRA) e análise de sequência genômica:

1. **PSMD6 (*Loriculus galgulus* e *Agapornis personatus*):**
   - O stop no aa 192 em *L. galgulus* é **polimórfico / heterozigoto** (21 reads com alelo funcional `CGA`/Arg e 20 reads com `TGA`/stop). O indivíduo preserva alelo íntegro.
   - O stop no aa 293 em *A. personatus* é corrigido por read independente (`SRR13483920`, AMNH) apresentando base `G` (`CGA`/Arg).
2. **ACD (*Loriculus*, *Psittacus erithacus*, *Psittacus timneh*, *Myiopsitta*):**
   - Em todas as espécies, o suposto stop prematuro cai exatamente na mesma coordenada do 1º éxon (aa 34–35).
   - O local situa-se em uma ilha CpG hiper-GC (**> 75 % de GC**), sofrendo drop-out severo de Illumina (cobertura 0× a 1×) e alinhamento imperfeito do miniprot. O restante do gene (éxons 2 a 15) possui ~30× de cobertura normal e integridade perfeita.
3. **Conclusão de Integridade:** **Nenhum dos genes do painel apresenta pseudogenização fixada em Psittaculidae**. A hipótese de perda gênica total é rejeitada pelos dados empíricos de sequenciamento.

---

## 5. Leitura Biológica Consolidada

1. **Rejeição da Perda Gênica Estrutural:** Não há deleção ou pseudogenização de freios tumorais/reparo em Psittaculidae.
2. **Divergência Seletiva na Maquinaria de Telômeros (M7):**
   - O sinal mais robusto e replicado em toda a árvore é o complexo de manutenção telomérica.
   - *TERT* opera sob seleção intensificada em Psittaculidae (\(k > 3\)), enquanto fatores associados (*DKC1* e *RTEL1*) mostram relaxamento seletivo.
   - No ramo terminal de *Melopsittacus*, praticamente todo o maquinário de *shelterin* e telomerase (*TERT*, *TERF1*, *TERF2*, *ACD*, *RTEL1*) sofreu surtos episódicos intensos de diversificação (\(q < 0,01\)), conversando diretamente com a alta taxa de neoplasias relatada em necropsias (20,8–23,2 %) em comparação com psitacídeos longevos.
3. **Modulação de Supressores Tumorais e Resposta ao Dano:**
   - Relaxamento pontual em *TSC2* (via mTOR) e no pró-apoptótico *BBC3* (PUMA), acompanhado de intensificação em antioxidantes centrais (*SOD1*) e resposta a quebras de fita dupla (*MRE11*).
