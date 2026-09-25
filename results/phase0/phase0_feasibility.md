# Relatório de viabilidade — Fase 0

**Data:** 25 set. 2026. **Reproduzir:** `bash scripts/00_fetch_external.sh && python scripts/01_build_phase0_tables.py && python scripts/02_phase0_analysis.py` (ver README).
**Status:** exploratório. Nenhuma análise molecular foi feita. Todos os contrastes abaixo derivam só do fenótipo e da árvore.

## 1. Cobertura fenotípica

| Item | Valor |
| --- | --- |
| Espécies na *supertree* (Burgio et al. 2019) | 413 (398 atuais + 15 extintas), altura 46,8 Ma |
| Espécies com expectativa de vida (Smeele et al. 2022, dados processados) | **218** |
| … com massa corporal e presentes na árvore | 218 (100 %) |
| Estimativas imprecisas (EP em log > 0,2) | 4 |
| Longevidade máxima compilada | **0**: AnAge/registros máximos não acessíveis deste ambiente |

**Desfecho principal proposto:** expectativa de vida **ao nascer**, em **cativeiro**, estimada por modelo Gompertz “bathtub” (BaSTA) sobre registros ZIMS, sexos agrupados. A cobertura (218 espécies, EP disponível para cada uma) é suficiente para a análise genômica, de modo que a longevidade máxima fica como **análise de sensibilidade** assim que for compilada, em coluna separada (`max_longevity_years`, hoje vazia).

Cautelas:
- Expectativa de vida ao nascer incorpora mortalidade juvenil e condições de cativeiro. Não mede taxa de senescência.
- Kakapo (*Strigops*) não tem estimativa porque não há registros zoológicos. Isso reduz Strigopidae a *Nestor* (2 espécies).

## 2. Alometria e sinal filogenético

`log(expectativa de vida) ~ log(massa)`, 218 espécies (`model_fits.tsv`):

| Modelo | Inclinação (EP) | Parâmetro | ΔAIC |
| --- | --- | --- | --- |
| OLS (sem filogenia) | 0,458 (0,025) | — | — |
| BM | 0,288 (0,043) | — | 69,9 |
| **λ de Pagel** | **0,349 (0,032)** | λ = 0,65 | **0** |
| OU | 0,382 (0,032) | α = 5,4 / altura da árvore | 14,9 |

- A inclinação é estável: excluir qualquer gênero com ≥ 5 espécies, excluir Cacatuidae ou excluir estimativas imprecisas mantém 0,34–0,38 (`sensitivity_allometry.tsv`).
- **Longevidade relativa** (resíduo do PGLS λ) tem sinal filogenético moderado: λ = 0,65, LR = 67,9, p ≈ 10⁻¹⁶. A expectativa de vida bruta tem λ = 0,92.

## 3. Onde está a longevidade relativa alta? (fenótipo, sem genes)

Gêneros com média de longevidade relativa > +0,25 em log (≈ +28 % acima do esperado pela massa), com ≥ 2 espécies:
- **Arini neotropicais**: *Psittacara* (+0,50), *Ara*, *Primolius*, *Anodorhynchus*, *Eupsittula*, *Pyrrhura*, *Aratinga*, e também *Pionites*;
- **Androglossini**: *Amazona* (+0,29; 21 espécies, amplitude −0,40 a +0,67);
- **Afro-malgaxes**: *Poicephalus*, *Coracopsis*;
- **Cacatuidae**: *Cacatua* (+0,30).

Esses grupos estão em regiões distantes da árvore: *Ara* × *Amazona* 25 Ma; Arini/*Amazona* × *Poicephalus* 32 Ma; Psittacidae × Cacatuidae 40,8 Ma (datas da *supertree*).

**Espécies prioritárias do protocolo** (`priority_species_phenotypes.tsv`), longevidade relativa (% acima/abaixo do esperado pela massa) e posição entre 218:

| Espécie | Expectativa de vida (anos, IC 95 %) | Relativa | Posição |
| --- | --- | --- | --- |
| *Ara macao* | 35,5 (33,0–38,2) | +86 % | 7 |
| *Cacatua galerita* | 25,2 (23,6–26,9) | +57 % | 25 |
| *Amazona aestiva* | 21,5 (20,3–22,8) | +55 % | 26 |
| *Anodorhynchus hyacinthinus* | 23,2 (22,1–24,4) | +10 % | 103 |
| *Cacatua moluccensis* | 18,2 (17,2–19,2) | +3 % | 116 |
| *Nymphicus hollandicus* | 8,3 (7,9–8,6) | 0 % | 127 |
| ***Nestor notabilis*** | 13,7 (12,7–14,7) | **−22 %** | 179 |
| *Melopsittacus undulatus* | 4,6 (4,5–4,7) | −25 % | 188 |

**Achado que muda o desenho:** com este desfecho, o **kea não é longevo para sua massa**, e a arara-azul-grande e a cacatua-das-molucas estão perto do esperado. Isso confirma a advertência do protocolo: rótulos “longevo” escolhidos por reputação não resistem ao ajuste por massa. Pode ser efeito do desfecho (mortalidade juvenil em cativeiro, manejo): a longevidade máxima precisa ser testada antes de concluir.

## 4. História evolutiva (reconstrução BM exploratória)

`ancestral_states_rel_logLE.tsv`, `largest_branch_changes.tsv`, `figures/fig1_tree_relative_LE.png`.

- Raiz: longevidade relativa ≈ 0 (EP 0,34), ou seja, compatível com um ancestral “mediano” e com um ancestral “longevo”.
- As **maiores mudanças estimadas estão em ramos terminais** (*Psittacula roseata*, *Cacatua sanguinea*, *Psittacara finschi*, *Amazona vittata*…). Isso é esperado sob BM: a reconstrução encolhe os nós internos em direção à média.
- **Nenhum ramo interno** tem mudança maior que ~1,5 EP. O maior aumento interno é no ancestral do clado Arini + *Pionites*/*Deroptyus* (47 espécies, +0,19, EP 0,13).
- Premissas não atendidas: λ = 0,65 indica que o BM puro é inadequado (AIC), e a reconstrução ignora o erro de medida das pontas.

**Leitura:** o padrão é compatível tanto com **vários aumentos independentes** (Arini, *Amazona*, *Poicephalus*, *Cacatua*) quanto com **ancestral relativamente longevo seguido de reduções** (Loriini, *Neophema*, *Cyanoramphus*, *Tanygnathus*, *Psittaculirostris*). A Fase 0 **não distingue** essas histórias. Modelos de mudança de regime (OU multi-ótimo, *l1ou*/*bayou*, RevBayes com erro de medida) são o próximo passo, com a árvore de Smith et al. 2023/2024 como sensibilidade.

## 5. Genomas × fenótipo

`assemblies.tsv`, `genome_phenotype_overlap.tsv`.

| Nível | Espécies com expectativa de vida |
| --- | --- |
| Montagem em escala cromossômica (VGP/GenomeArk, verificada no S3) | **7**: *Amazona ochrocephala*, *Ara ararauna*, *Cacatua galerita*\*, *Guaruba guarouba*, *Lathamus discolor*, *Melopsittacus undulatus*, *Psittacula krameri* |
| Draft ou VGP em andamento | 7: *Amazona aestiva*, *A. vittata*, *Ara macao*, *Anodorhynchus hyacinthinus*, *Nestor notabilis*, *Agapornis roseicollis*, *Psitteuteles goldiei* |
| Accession relatada, nível desconhecido | 3: *Cacatua sanguinea*, *C. tenuirostris*, *C. ducorpsii* |
| **Total mapeado** | **17 de 218** |

\* *C. galerita*: accession GCA_035583095.1 citada em resumo de busca; não verificada no NCBI.

**Limitação central:** o NCBI Datasets estava bloqueado. O lote de **Hains et al. 2022 (94 espécies)** e o de 22 espécies (F1000 2020) **ainda não foram mapeados espécie a espécie**, então 17 é um **piso**. A sobreposição real pode ser várias vezes maior, mas só o mapeamento dirá quanto. Métricas BUSCO/N50 só existem hoje para 4 montagens VGP (N50 de scaffold de 83–110 Mb).

## 6. Poder (simulação)

`power_simulation.tsv`: associação entre dois caracteres coevoluindo por BM com correlação evolutiva *r*, testada por PGLS; 400 simulações por célula. É um **proxy** otimista para um teste gene a gene (taxas gênicas têm mais ruído que um caráter).

| Amostra | r = 0,2 | r = 0,4 | r = 0,6 |
| --- | --- | --- | --- |
| 17 spp. atuais, α = 0,05 | 0,13 | 0,40 | 0,77 |
| 17 spp., α = 5×10⁻⁶ (≈ Bonferroni para 10⁴ genes) | 0,00 | 0,00 | 0,01 |
| 40 spp. aleatórias, α = 5×10⁻⁶ | 0,00 | 0,01 | 0,35 |
| 94 spp. aleatórias, α = 5×10⁻⁶ | 0,01 | 0,30 | 0,99 |
| 218 spp., α = 5×10⁻⁶ | 0,05 | 0,94 | 1,00 |

- Com 17 espécies, uma **varredura genômica ampla não tem poder**. Mesmo o painel dirigido (~100 genes) só detectaria efeitos grandes.
- Com ~94 espécies, efeitos moderados a grandes passam a ser detectáveis em escala genômica.

## 7. Eixos integrados — viabilidade resumida

- **Câncer (6A):** o contraste de *Melopsittacus* existe na literatura [C0, C1], mas os comparadores publicados são “outros psitacídeos” ou “outras aves”. Os **parentes reais** do periquito são os **Loriini** (~9,75 Ma; 61 espécies a ≤ 20 Ma, várias com expectativa de vida e *Psitteuteles goldiei* com genoma draft). A calopsita (41 Ma) é controle de manejo, não parente próximo. Não há dados de necropsia com denominador para lóris neste levantamento.
- **Neurodegeneração (6B):** nenhum estudo de neuropatologia do envelhecimento espontâneo foi localizado em nenhuma espécie. O que existe é doença genética precoce [N1] e plasticidade vocal em adultos jovens [N2]. O entregável viável agora é a **matriz de lacunas**, não uma análise.
- **Oxidação (6C):** evidência heterogênea mas utilizável: produção de ROS semelhante [R2]; GPx maior, dano não menor [R3]; resistência celular induzida em periquito [R1]; telômeros longos com atrito mais rápido em longevos (longitudinal) [L13]. Isso já rejeita a versão simples “mais antioxidante = menos dano”.

## 8. Decisão

> **Cenário B agora, com caminho definido para A/B em escala.**

- **B:** o fenótipo contínuo é informativo (218 espécies, sinal filogenético claro, alometria robusta), mas a história ancestral é incerta e não há transições internas bem sustentadas. **Priorizar associação genótipo–fenótipo contínuo; não afirmar múltiplas origens.**
- **Não C para o fenótipo, mas C para a genômica com a amostra atual:** 17 espécies mapeadas dão poder insuficiente para varredura genômica. O painel VGP (7–9 espécies) serve como **piloto técnico** (extração, ortologia, alinhamento, TERT).
- **Condição para escalar:** mapear o lote Hains et al. (94 spp.) + F1000 (22 spp.) ao fenótipo. Se a sobreposição ficar ≥ ~60 espécies com BUSCO aceitável e distribuída pelos clados longevos e não longevos, seguir com RERconverge em escala.
- **Condição para testar convergência (A):** um modelo de regimes com erro de medida precisa identificar ≥ 3 mudanças para “longevidade relativa alta” em clados independentes, com suporte robusto às duas árvores.

## 9. Bloqueios deste ambiente

Hosts negados pela política de rede: NCBI (Datasets, E-utilities, FTP), Crossref, Europe PMC, OpenAlex, Dryad, Zenodo, Figshare, bioRxiv, Ensembl, BUSCO, CRAN. GitHub e GenomeArk S3 funcionaram. Para completar a Fase 0, rodar `00_fetch_external.sh` e a auditoria NCBI em um ambiente com esses domínios liberados.
