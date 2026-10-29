# Rascunho — Resultados (28 out. 2026)

Rascunho para o texto do trabalho. Números conferidos nos relatórios do repositório (`results/phase1_cancer/cancer_report.md`, `results/phase1_audit/audit_report.md`, `results/phase2_pilot*/`). A análise principal de genes (58 espécies) ainda está rodando; as seções 3 e 4 trazem o piloto, que é exploratório.

## 1. Neoplasias: o periquito é um caso extremo entre os papagaios

Em necropsias de zoológicos credenciados pela AZA (Kapsetaki et al. 2024; 5.729 necropsias de 108 espécies de aves), o periquito-australiano (*Melopsittacus undulatus*) teve neoplasia em 20,8 % (99/477) e malignidade em 17,4 %. Nos lóris (Loriini, os parentes mais próximos com dados), os valores foram 1,5 % (6/399) e 0,8 % (razão de chances 17; Fisher, p ≈ 10⁻²¹). Considerando a variação entre as 15 espécies de papagaios, o periquito ficou fora da distribuição das outras 14 (p ≈ 0,002, modelo beta-binomial). A calopsita (*Nymphicus hollandicus*), também domesticada, teve 11,4 % (8/70).

O padrão se repetiu numa base independente, de clínicas e criadores: Drury R. Reavill Pathology Database, UC Davis, 10.152 necropsias de papagaios com ≥ 3 tecidos. Ali o periquito teve tumor em 23,2 % (n = 652; 25,9 % padronizado por idade) e os lóris em 5,0 % (n = 399; 5,4 %), razão de ~4,6× bruta e ~4,8× padronizada (p ≈ 10⁻¹⁶). Araras, papagaio-cinzento, cacatuas e conures ficaram entre 3 % e 6 %. A calopsita (9,1 %, n = 1.348) ficou bem abaixo do periquito (p ≈ 10⁻¹⁷), o que enfraquece a domesticação como explicação única. A diferença entre bases (14× contra 5×) pode refletir população (zoológico × pets), idade à morte e critérios diagnósticos.

## 2. Auditoria dos genomas: controle de amostras falho no lote Iridian

Das 316 montagens de Psittaciformes no NCBI, 235 são da Iridian Genomes; 96 cobrem espécies com dados de longevidade. Todas foram auditadas quanto a completude, identidade do espécime e contaminação.

- **Completude** (compleasm, aves_odb10): mediana de 76 % de genes BUSCO completos. 57 de 96 passam no critério revisto (≥ 70 %); só 10 passariam no critério original (≥ 85 % e N50 de contig ≥ 20 kb). Nenhuma montagem com N50 < 5 kb (31) passou.
- **Identidade:** 61 das 96 montagens contêm o mitogenoma completo de outra espécie, em 88 % dos casos uma espécie também sequenciada pela Iridian, o que indica contaminação cruzada entre amostras. Pelo genoma nuclear (distâncias Mash), quatro montagens não são a espécie declarada: "*Agapornis nigrigenis*" é um passeriforme; "*Cacatua haematuropygia*" e "*Trichoglossus haematodus*" estão misturadas com esse mesmo genoma; e "*Aratinga jandaya*" é um *Poicephalus meyeri*.
- **Conjunto final:** com as regras registradas no plano antes de qualquer análise de genes (um genoma por espécie, BUSCO ≥ 70 %, identidade não reprovada, preferência por montagens cromossômicas de outras fontes), a lista principal ficou com 58 espécies, 15 de Psittaculidae.

## 3. Piloto da H8 (exploratório): 101 genes, 17 genomas de alta qualidade

Foram testados 94 genes dos oito módulos (M1–M8) em 17 espécies com montagens de leitura longa ou cromossômicas: 10 Psittaculidae e 7 de outras famílias. Nenhum lóri tem montagem desse tipo.

- **Sem padrão modular.** A intensidade de seleção (RELAX, *k*) não diferiu entre módulos (permutação, p = 0,67). Não houve relaxamento generalizado em Psittaculidae, nem perda gênica restrita à família. As interrupções de leitura encontradas aparecem dentro e fora da família e se concentram nas montagens mais sujeitas a erro.
- **Genes isolados em Psittaculidae** (BH, q < 0,10): relaxamento em CASP8 (apoptose; *k* ≈ 0), DKC1, TXN e NFE2L2; seleção intensificada em **TERT** (*k* = 3,2; q = 8 × 10⁻⁶), POT1, ACD, SQSTM1, HSF1, KEAP1, SOD1, SOD2 e CAT.
- **Ramo do periquito** (aBSREL): seleção episódica em **TERT** (p = 2 × 10⁻⁵), **TERF2** (p = 3 × 10⁻⁴) e **RTEL1** (p = 3 × 10⁻¹⁴). No tronco de Psittaculidae, ACD (p = 5 × 10⁻⁵). Os três genes do ramo do periquito já apareciam no piloto menor, de 8 espécies.

## 4. As mudanças nos genes de telômero do periquito se repetem em outro indivíduo

Nos genes de telômero (M7), 5,9 % dos aminoácidos são exclusivos do periquito, contra 0,7 % nos demais módulos (mediana; análise feita após ver os dados). Para descartar erro de montagem, a análise foi repetida num segundo periquito, sequenciado por outra tecnologia (Iridian, Illumina). Ele carrega 61 das 63 substituições exclusivas do RTEL1 (97 %) e 25 das 26 do POT1 (96 %). Por outro lado, a interrupção do MRE11 no genoma de referência (VGP) não aparece no segundo indivíduo, então era erro de montagem. TP53 e TERT não puderam ser comparados, porque esses genes vieram fragmentados na montagem Iridian.

## 5. Leitura e limites

- Pela tabela de interpretação pré-registrada, a ausência de padrão modular e de relaxamento generalizado **enfraquece a versão "perda gênica" da H8**. As versões funcionais (H8a/H8b) seguem em aberto. Sinais concentrados no ramo do periquito, e não no tronco de Psittaculidae, **enfraquecem a H8c**.
- O sinal mais consistente é a **evolução acelerada da maquinaria de telômeros no periquito**, a espécie com mais neoplasias. Ela foi confirmada em dois indivíduos e por teste de seleção, e é compatível com o papel da telomerase em câncer. Não se sabe o efeito funcional dessas mudanças.
- Limites: o piloto tem 17 espécies e nenhum lóri, então não separa o periquito do ancestral periquito + lóris. As montagens do piloto ampliado não passaram pela auditoria BUSCO. A métrica de "mudanças exclusivas" foi escolhida depois de ver os dados. A análise principal (58 espécies) e a extensão pré-declarada são os testes confirmatórios.
