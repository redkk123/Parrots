# Plano de análise (fixado antes da análise molecular)

**Versão:** 1.0, 25 set. 2026. Alterações posteriores devem ser registradas abaixo, em “Histórico”, com data e justificativa, **antes** de rodar os testes que elas afetam.

## 1. Desfecho

- **Principal:** `log_life_expectancy`, a expectativa de vida ao nascer em cativeiro (Smeele et al. 2022), com EP por espécie.
- **Covariável obrigatória:** `log_body_mass`.
- **Forma do fenótipo nas análises moleculares:** o caráter contínuo `log_life_expectancy` com `log_body_mass` como covariável, **ou** a longevidade relativa (resíduo do PGLS λ, inclinação 0,349). Se usar resíduos, propagar a incerteza por reamostragem dos coeficientes do PGLS e do EP das pontas.
- **Sensibilidade 1:** longevidade máxima documentada (coluna separada, a compilar), com o mesmo ajuste por massa.
- **Sensibilidade 2:** excluir espécies com EP(log) > 0,2 (4 espécies).
- **Não usar** “longevo” binário como desfecho principal. Categorias ficam para análises secundárias.

## 2. Covariáveis e justificativa causal

| Variável | Papel | Decisão |
| --- | --- | --- |
| Massa corporal | confundidor (afeta longevidade e taxas moleculares) | incluir sempre |
| Tamanho cerebral relativo | possível mediador ou causa comum | **não** ajustar no modelo principal; análise secundária |
| Tamanho da postura, idade de 1ª reprodução, desenvolvimento | história de vida (mediador/confundidor) | só em sensibilidade, dada a cobertura incompleta |
| Qualidade da montagem (BUSCO, N50 de contig) | confundidor técnico | covariável ou filtro; testar associação com o fenótipo antes |

## 3. Árvores

- **Principal:** *supertree* de Burgio et al. 2019 (a mesma de Smeele et al. 2022).
- **Sensibilidade obrigatória:** árvore datada de Smith et al. 2023/2024 (UCE), após harmonização taxonômica.
- Para taxas gênicas (RERconverge), comprimentos de ramo estimados por gene sobre a topologia de espécies fixada.

## 4. Critérios de inclusão de genomas (fixar antes dos testes)

- BUSCO (aves_odb10) completo ≥ 85 % **e** N50 de contig ≥ 20 kb para o conjunto amplo; ≥ 95 % e escala cromossômica para o painel piloto.
- Um genoma por espécie. Preferir RefSeq ou VGP curado.
- Ausência ou truncamento de gene não conta como perda sem reanotação independente (miniprot/TOGA).
- Testar se BUSCO e N50 se associam à longevidade relativa (PGLS). Se sim, incluí-los como covariáveis.

## 5. Teste principal e família de hipóteses

1. **Principal:** RERconverge, correlação entre taxas relativas por gene e o fenótipo contínuo (com massa), permulações filogenéticas para calibrar p; FDR (BH) em todos os genes elegíveis.
2. **Painel dirigido (família separada):** os ~100 genes da tabela do protocolo, com FDR própria. Declarado como teste de hipótese prévia, não enriquecimento.
3. **Seleção (complementar):** aBSREL/BUSTED nos ramos com longevidade relativa alta **definidos pelo fenótipo** (item 7); RELAX para separar relaxamento de intensificação. Concordância entre ferramentas no mesmo alinhamento não conta como replicação.
4. **Vias:** enriquecimento com background = genes elegíveis testados; controle por tamanho e sobreposição de vias; permutações que preservem a filogenia.

## 6. Robustez obrigatória

- Leave-one-clade-out: Arini, *Amazona*, Cacatuidae, Loriini, Psittaculini.
- Sensibilidade à árvore (item 3) e ao filtro de alinhamento (antes/depois).
- Reportar tamanhos de efeito com IC, não só p.

## 7. Contrastes e ramos de interesse (definidos agora, só pelo fenótipo)

Ramos “longevidade relativa alta” para testes por ramo, **exploratórios**:
- ancestral de Arini + *Pionites*/*Deroptyus* (maior aumento interno estimado: +0,19, EP 0,13);
- ramos terminais das espécies com genoma e longevidade relativa ≥ +0,4: *Ara macao*, *Amazona vittata*, *A. ochrocephala*, *A. aestiva*, *Cacatua galerita*, *C. sanguinea* (se o genoma for confirmado).

Contrastes próximos (pares irmãos ou quase, com genoma nos dois lados, a completar com o lote Hains):
- *Ara macao* (+0,62) × *Ara ararauna* (+0,23);
- *Amazona ochrocephala*/*aestiva* × *Amazona* de longevidade relativa baixa (a identificar entre as espécies com genoma);
- *Melopsittacus* (−0,29) × Loriini (*Psitteuteles goldiei*, −0,16; outros a mapear).

Esses rótulos **não** implicam origens independentes. A convergência só será testada se o modelo de regimes (feasibility §8) sustentar ≥ 3 mudanças independentes.

## 8. Eixos integrados

- **Câncer:** comparar *Melopsittacus* primeiro com Loriini, depois com Psittaculidae (~32 Ma), e usar *Nymphicus* como controle de manejo. Modelo binomial/beta-binomial hierárquico (espécie, fonte) só com denominadores. Separar benigno/maligno e tecido.
- **Neuro:** entregável = matriz de lacunas; nenhuma inferência de proteção.
- **Oxidação:** tabelar por etapa (produção, neutralização, vulnerabilidade, reparo, tolerância). Não combinar ensaios incompatíveis em um índice.

## 9. Recursos

Piloto (7–9 genomas VGP, ~100 genes): CPU, ≤ 12 GB RAM. Medir tempo e memória antes de escalar para ~100 genomas.

## Histórico

- 2026-09-25: v1.0.
- 2026-09-25 (Fase 0b, antes de qualquer análise molecular):
  - Amostra genômica: camadas A–C (105 spp. com expectativa de vida) em [`../phase0b/species_data_matrix.tsv`](../phase0b/species_data_matrix.tsv). Camada D (5–20×) só em sensibilidade.
  - A convergência (H3) **não** será testada com os rótulos do §7: a análise de regimes não sustenta aumentos independentes (ver `phase0b_report.md` §7).
  - Adicionado teste dirigido no ramo-tronco de Psittaculidae (aBSREL/RELAX), com a ressalva de evento único.
  - Adicionadas análises contínuas dentro de Psittaculidae e dentro do restante.
  - Critério de robustez: mesma direção com a longevidade máxima ajustada pelo tamanho amostral (AnAge). Espécies com `qc_outcomes_inconsistent` só em sensibilidade.
  - Contrastes irmãos: usar [`../phase0b/sister_pair_contrasts.tsv`](../phase0b/sister_pair_contrasts.tsv). Pares em que os desfechos discordam são marcados como ambíguos.
- 2026-09-26 — **H8 (travas moleculares modulares, retenção parcial em Psittaculidae)** adicionada ao protocolo (`docs/protocol.md`), antes de qualquer análise molecular. Módulos e testes fixados abaixo; alterações exigem nova entrada datada.

## 10. Teste de H8 (fixado em 2026-09-26)

### Módulos (listas congeladas; presença e ortologia a confirmar em aves antes dos testes)

| Módulo | Função | Genes |
| --- | --- | --- |
| M1 | resposta a dano / checkpoints | ATM, ATR, CHEK1, CHEK2, TP53, MDM2, MDM4, CDKN1A, BUB1B, BUB3, MAD2L1 |
| M2 | controle proliferativo | RB1, CDKN2A/CDKN2B (locus a confirmar em aves), CDKN1B, PTEN, CCNE1, MYC, TSC1, TSC2 |
| M3 | apoptose / senescência | BAX, BAK1, BCL2, BCL2L1, MCL1, BBC3, PMAIP1, APAF1, CASP3, CASP8, CASP9 |
| M4 | reparo de DNA (inclui lesões oxidativas) | BRCA1, BRCA2, RAD51, MRE11, NBN, RAD50, PARP1, ERCC3, POLK, OGG1, MUTYH, NUDT1, APEX1, POLB, XRCC1 |
| M5 | resposta **induzida** ao estresse | NFE2L2, KEAP1, HSF1, HMOX1, NQO1, SQSTM1, TXNRD1, GCLC, GCLM, SRXN1 |
| M6 | antioxidantes basais | SOD1, SOD2, SOD3, CAT, GPX1, GPX4, PRDX1–PRDX6, TXN, GSR |
| M7 | telômeros | TERT, POT1, TERF1, TERF2, TINF2, ACD, TERF2IP, RTEL1, DKC1 (TERC em análise própria) |
| M8 | proteostase / autofagia | PSMD6, subunidades PSMA/PSMB, ATG5, ATG7, BECN1, PINK1, PRKN, BNIP3, BNIP3L, TFEB |

### Testes

1. **Por gene:** RELAX (intensificação × relaxamento) com ramos de teste = tronco de Psittaculidae + ramos internos de Psittaculidae e referência = restante; aBSREL no tronco. Perda/pseudogenização só com reanotação independente (miniprot/TOGA contra a referência cromossômica de *Melopsittacus* e de uma espécie longeva) e cobertura adequada.
2. **Por módulo:** proporção de genes com relaxamento significativo (FDR 10 % dentro do painel) e média do parâmetro *k* do RELAX. Nulo: 1.000 conjuntos aleatórios de genes do universo elegível, pareados por tamanho de alinhamento e dN/dS de fundo.
3. **Previsão central (1–2):** heterogeneidade entre módulos, com M1–M3 mais relaxados que M5. Teste de contraste pré-especificado: média(*k*) de M1–M3 < média(*k*) de M5, por permutação dos rótulos de módulo. Se todos os módulos estiverem igualmente relaxados, H8 perde a parte “parcial”.
4. **Previsão 3:** RERconverge com fenótipo contínuo **dentro** de Psittaculidae e **dentro** do restante. Um sinal que só aparece entre clados é tratado como evento único.
5. **Previsão 4:** só se houver ≥ 15 espécies com necropsias e denominadores: prevalência relativa de neoplasia ~ número de módulos sem relaxamento + massa (modelo binomial filogenético).
6. **Controles:** qualidade da montagem/cobertura como covariável; *Melopsittacus* e Loriini testados também como ramos terminais separados; repetição com a árvore de Smith et al.

### Interpretação pré-registrada

- **Apoio a H8:** heterogeneidade entre módulos na direção prevista, localizada em Psittaculidae, robusta à qualidade dos dados e à árvore.
- **Contra H8:** nenhum módulo diferente, relaxamento uniforme, ou módulos alterados sem relação com controle proliferativo.
- **Inconclusivo:** poucos genes testáveis por módulo (< 6 após controle de qualidade) ou efeitos na direção prevista sem significância.
