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
- 2026-09-26 (revisão da H8 após crítica externa, antes de qualquer análise molecular): H8 dividida em **H8a** (proteção modular), **H8b** (efeito duplo: células danificadas sobrevivem mais) e **H8c** (história: mudanças no tronco de Psittaculidae). Retirada a previsão “M5 conservado / controle proliferativo perdido” como consequência de R1: NRF2 e resposta ao choque passam a ser candidatos. Genes intactos deixam de refutar a hipótese funcional e passam a enfraquecer só a versão de perda gênica. A distribuição dos tumores por tecido deixa de contar como evidência. A seção 10 foi reescrita.

## 10. Teste de H8 (fixado em 2026-09-26, revisado no mesmo dia)

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

Nenhum módulo tem direção prevista a partir de R1. Sob H8a, espera-se diferença de eficácia entre módulos de **prevenção/reparo** (M4–M6, M8) e de **controle de células danificadas** (M1–M3). Sob H8b, espera-se, em particular, limiar mais alto em M3 (e possivelmente M1). Qual módulo difere é uma pergunta, não uma premissa.

### Nível 1 — genômico (este projeto, in silico)

1. **Perda gênica:** presença, integridade da ORF e pseudogenização por reanotação independente (miniprot/TOGA contra *Melopsittacus* e uma espécie longeva com montagem cromossômica), com cobertura de reads verificada. Só esse resultado testa a **versão perda gênica**.
2. **Restrição seletiva:** RELAX (tronco de Psittaculidae + ramos internos × restante), aBSREL no tronco, RERconverge contínuo **dentro** de Psittaculidae e **dentro** do restante. Resultado por módulo: proporção de genes com relaxamento (FDR 10 % dentro do painel) e média de *k*. Nulo: 1.000 conjuntos aleatórios do universo elegível, pareados por tamanho de alinhamento e dN/dS de fundo.
3. **Heterogeneidade entre módulos:** teste bicaudal por permutação dos rótulos de módulo. Heterogeneidade é compatível com modularidade (H8a ou H8b); ausência de heterogeneidade não refuta a versão funcional.
4. **Localização (H8c):** sinais concentrados no tronco de Psittaculidae × espalhados pela árvore × restritos a *Melopsittacus*. Por ser evento único, sinal no tronco não é atribuído à longevidade sem suporte dentro dos clados.

### Nível 2 — regulatório (se houver dados comparáveis; senão, relatar a lacuna)

Expressão basal e **induzida** por estresse dos genes dos módulos em tecidos/linhagens comparáveis, elementos não codificantes conservados próximos aos genes, e variantes de splicing. Testa a versão **menor eficácia regulatória**.

### Nível 3 — celular (proposta para validação futura)

Experimento da tabela H8a × H8b no protocolo: sobrevivência, dano inicial e **residual**, parada do ciclo em células com dano, senescência/apoptose para dano equivalente e destino das sobreviventes. Espécies: *Melopsittacus*, 1–2 Loriini, 1 espécie longeva fora de Psittaculidae, codorna. É o único nível que distingue H8a de H8b.

### Interpretação pré-registrada

| Resultado | Leitura |
| --- | --- |
| Perda gênica confirmada em módulos de controle, localizada em Psittaculidae | apoia H8 (versão perda gênica) e H8c |
| Genes intactos, sem relaxamento | enfraquece a versão perda gênica; H8a/H8b funcionais ficam **em aberto** (níveis 2–3) |
| Relaxamento heterogêneo entre módulos | compatível com modularidade; não indica direção de eficácia nem distingue H8a/H8b |
| Sinais fora de Psittaculidae ou só em *Melopsittacus* | enfraquece H8c, sem afetar H8a/H8b |
| Ensaios celulares sem diferença entre periquito/Loriini e espécies longevas | refuta a versão funcional |
| Dano residual maior + menor parada/eliminação para o mesmo dano | favorece H8b |
| Dano residual menor + controle de células danificadas preservado ou reduzido | favorece H8a |

Tumores por tecido (renais, pituitários, proventriculares, lipomatosos) servem só para gerar candidatos, não como evidência. *Macrorhabdus*, hormônios, dieta e obesidade são explicações concorrentes específicas por tecido.
