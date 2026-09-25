# Fase 0b — Mapeamento genoma × fenótipo, desfecho alternativo e história da longevidade

**Data:** 25 set. 2026. **Continua:** [`../phase0/phase0_feasibility.md`](../phase0/phase0_feasibility.md).
**Reproduzir:** `scripts/03_sra_parrot_inventory.py` → `04_build_data_matrix.py` → `05_phase0b_analysis.py`, e depois `scripts/R/validate_phase0.R` e `scripts/R/regime_shifts.R` (ambientes: `scripts/setup_r_envs.sh`).

A Fase 0 terminou com três pendências: (1) mapear os genomas short-read (Hains et al. 2022 + outros) contra o fenótipo, (2) compilar a longevidade máxima como desfecho de sensibilidade e (3) testar se a história da longevidade relativa sustenta mudanças independentes. Esta fase resolve as três na medida em que os dados públicos permitem.

## 1. Como os dados foram obtidos sem NCBI Datasets

O NCBI Datasets e o E-utilities continuam bloqueados neste ambiente. Porém o NCBI publica **todos os metadados do SRA** como Parquet no programa AWS Open Data (`s3://sra-pub-metadata-us-east-1`, snapshot de 25 set. 2026), e esse host é permitido. A coluna `organism` dos 30 arquivos (13,6 GB) foi varrida por leitura de intervalos HTTP, o que custou ~1,1 GB de download, e as demais colunas só foram lidas nos grupos de linhas com papagaios.

- 3.300 corridas de Psittaciformes; 2.052 rotuladas “WGS”.
- **Excluídas** as corridas “WGS” com seleção de biblioteca *Hybrid Selection* (543, sobretudo captura de UCEs do AMNH/LSU) e *Restriction Digest* (1). Captura dirigida não é cobertura genômica.
- Nomes resolvidos para a árvore em camadas: binômio (trinômios colapsados), sinônimos do banco de Burgio e tabela curada ([`taxonomy_synonyms.tsv`](../../data/curated/taxonomy_synonyms.tsv)), que distingue **mudança de gênero** (mesma espécie) de **espécie separada depois** (aproximação; não usada para a camada principal). Híbridos e “sp.” foram descartados.
- Cobertura estimada = bases do melhor indivíduo (BioSample) / 1,2 Gb. É estimativa, não cobertura medida por mapeamento.

## 2. Camadas de dados genômicos × fenótipo

[`species_data_matrix.tsv`](species_data_matrix.tsv): uma linha por espécie da árvore (413).

| Camada | Critério | Com expectativa de vida | Sem |
| --- | --- | --- | --- |
| A | montagem cromossômica VGP (GenomeArk) | 6 | 2 |
| B | long-read (PacBio/ONT) ≥ 15× num indivíduo | 14 | 2 |
| C | short-read ≥ 20× num indivíduo | 85 | 34 |
| D | short-read 5–20× | 39 | 32 |
| E | short-read < 5× | 4 | 1 |
| — | sem WGS público | 70 | 124 |

- **Utilizáveis para análise gene a gene (A–C): 105 espécies com expectativa de vida** (93 Psittacidae, 11 Cacatuidae, 1 Strigopidae) e 88 com longevidade máxima. Na Fase 0 eram 17.
- **96 das 105 vêm do lote de Hains et al. 2022** (Iridian Genomes, um BioProject por espécie, um indivíduo cada). É a base do projeto, e sua qualidade (BUSCO, contaminação, identificação do espécime) precisa ser auditada antes de qualquer teste.
- Das 55 espécies no quartil superior de longevidade relativa, 34 são utilizáveis. A amostra utilizável cobre quase toda a amplitude do fenótipo (DP 0,31 vs. 0,35 no total).

## 3. Poder com o conjunto real de espécies

[`power_realised_sets.tsv`](power_realised_sets.tsv): mesma simulação da Fase 0 (caracteres coevoluindo por BM, PGLS), agora na árvore das espécies realmente disponíveis.

| Conjunto | n | r = 0,3 | r = 0,4 | r = 0,6 |
| --- | --- | --- | --- | --- |
| A–C, α = 0,05 | 105 | 0,88 | 0,99 | 1,00 |
| A–C, α = 5×10⁻⁶ (genoma inteiro) | 105 | 0,08 | 0,37 | 0,99 |
| A–D, α = 5×10⁻⁶ | 144 | 0,16 | 0,61 | 1,00 |

Uma varredura genômica agora é viável para efeitos grandes. Para efeitos moderados, o painel dirigido e os testes por via continuam sendo o caminho de maior poder.

## 4. O desfecho importa muito: expectativa de vida × longevidade máxima

Longevidade máxima: AnAge (com origem, classe de tamanho amostral e qualidade; registros de qualidade “low” excluídos) e, para lacunas, o Amniote Life History Database. Os dois bancos **não são independentes**. [`max_longevity_sources.tsv`](max_longevity_sources.tsv).

- Alometria semelhante: inclinação 0,32 (EP 0,03), λ = 0,50, 171 espécies.
- **O recorde cresce com o tamanho amostral:** +0,13 em log por classe (EP 0,03). É o esperado para um valor extremo: quem tem milhões de indivíduos em cativeiro (periquito, calopsita) acumula recordes maiores. Por isso a longevidade máxima relativa foi recalculada **ajustando pela classe de tamanho amostral** (152 espécies).
- **Concordância fraca entre as medidas:** Spearman ρ = 0,39 (0,35 após ajuste; 138 espécies). Em 15 pares irmãos com as duas medidas, só **7** concordam sobre qual espécie é a mais longeva.

[`priority_species_two_outcomes.tsv`](priority_species_two_outcomes.tsv), [`figures/fig3_two_outcomes.png`](figures/fig3_two_outcomes.png):

| Espécie | Expectativa de vida relativa (posição/218) | Longevidade máxima relativa, ajustada (posição/152) |
| --- | --- | --- |
| *Nymphicus hollandicus* | 0 % (127) | +46 % (9) |
| *Cacatua moluccensis* | +3 % (116) | +38 % (16) |
| *Cacatua galerita* | +57 % (25) | +31 % (21) |
| *Amazona aestiva* | +55 % (26) | +28 % (24) |
| *Psittacus erithacus* | +28 % (73) | +24 % (27) |
| *Melopsittacus undulatus* | −25 % (188) | +16 % (41) |
| *Nestor notabilis* | −22 % (179) | +13 % (44) |
| *Anodorhynchus hyacinthinus* | +10 % (103) | +2 % (56) |
| *Ara macao* | +86 % (7) | −36 % (129) |

**Consequência para a Fase 0:** a conclusão “o kea não é longevo para sua massa” **depende do desfecho**. Pela longevidade máxima ajustada, kea, periquito e calopsita ficam acima do esperado.

**Controle de qualidade:** em 15 espécies o recorde é < 1,3 × a expectativa de vida média, e na *Ara macao* ele é até **menor** que a média (33 vs. 35,5 anos). Isso é internamente impossível e aponta erro em ao menos uma das fontes. As espécies afetadas são sobretudo **araras e amazonas**, exatamente os clados que a expectativa de vida coloca no topo. Há duas hipóteses concorrentes: (a) recordes do AnAge desatualizados; (b) extrapolação do modelo de sobrevivência (Gompertz “bathtub” sobre dados muito censurados) inflando a expectativa de vida das espécies mais longevas. Coluna `qc_outcomes_inconsistent` na matriz.

## 5. Contrastes irmãos com genoma dos dois lados

[`sister_pair_contrasts.tsv`](sister_pair_contrasts.tsv): 29 pares irmãos na árvore podada às 105 espécies utilizáveis. Por construção, os pares são filogeneticamente independentes entre si. Oito têm diferença ≥ 0,3 em log de expectativa de vida relativa, por exemplo *Cacatua sanguinea* × *C. haematuropygia* (11,6 Ma; Δ = 0,77, e a longevidade máxima concorda) e *Poicephalus senegalus* × *P. gulielmi* (21,8 Ma; concorda). Pares em que as duas medidas discordam (*Cacatua alba* × *C. moluccensis*, *Ara macao* × *A. chloropterus*) devem ser tratados como ambíguos.

## 6. Validação das análises em Python contra R

[`validation_python_vs_R.json`](validation_python_vs_R.json). O R 4.4 foi instalado via conda-forge.

- PGLS λ: Python, `nlme::gls(corPagel)` e `phylolm` dão λ = 0,6534 e inclinação 0,3491 idênticos. O EP da inclinação passou a usar a variância residual não enviesada (0,0322, como `nlme`).
- AIC BM / λ / OU idênticos ao `phylolm` (172,0 / 102,1 / 117,0), e α do OU idêntico.
- Estados ancestrais: `phytools::fastAnc` × Python nos 191 nós internos: diferença máxima 8×10⁻⁵, correlação das estimativas e dos EP = 1.

## 7. Mudanças de regime na longevidade relativa (PhylogeneticEM)

[`regime_shifts.json`](regime_shifts.json), [`figures/fig4_regime_shifts_rel_LE.png`](figures/fig4_regime_shifts_rel_LE.png), [`shift_clade_contrast.tsv`](shift_clade_contrast.tsv).

O modelo usado é um OU com ótimos que mudam em ramos desconhecidos (scOU), com busca por EM de 0 a 12 mudanças e seleção de *K* por LINselect (`BGHuni`). Os critérios DDSE e Djump não estão disponíveis para um único caráter nesta versão (1.8.1). O erro de medida das pontas não é modelado.

**Expectativa de vida relativa (218 espécies): K = 3, todas reduções.**

| Ramo | Mudança no ótimo (log) | Espécies no ajuste |
| --- | --- | --- |
| ancestral de **Psittaculidae** (lóris, *Psittacula*, *Agapornis*, *Loriculus*, rosellas, *Neophema*, *Cyclopsitta*, *Eclectus*, *Melopsittacus*…) | **−0,35** | 101 |
| *Psittaculirostris desmarestii* (ponta) | −1,30 | 1 |
| *Tanygnathus megalorynchos* (ponta) | −1,36 | 1 |

- Ótimo da raiz: **+0,22**. O regime ancestral é “relativamente longevo” e continua em Cacatuidae, Strigopidae (*Nestor*) e nos Psittacidae neotropicais e afro-tropicais.
- As duas mudanças de ponta isoladas são mais parcimoniosamente **erro de estimativa** do que evolução.
- α = 25,6 por altura da árvore (meia-vida ≈ 1,3 Ma): as espécies ficam perto do ótimo do seu regime, e o modelo funciona quase como um agrupamento de regimes.

**Longevidade máxima ajustada pela amostra (152 espécies): K = 0.** O critério não sustenta nenhuma mudança.

**Teste direto do clado (PGLS λ com indicador de Psittaculidae):**

| Desfecho | n | Efeito | p |
| --- | --- | --- | --- |
| Expectativa de vida relativa | 218 | −31 % | 0,003 |
| Longevidade máxima relativa | 171 | −23 % | 0,04 |
| Longevidade máxima ajustada pela amostra | 152 | −18 % | 0,03 |

Retirar cada subclado (Loriini, Platycercini + Pezoporini, Psittaculini, *Agapornis* + *Loriculus*, Cyclopsittini, *Melopsittacus*) mantém o efeito entre −28 % e −35 % (p ≤ 0,015). **É o único padrão que sobrevive à troca de desfecho.**

**Como ler isso:**
- A história mais bem sustentada pelos dados atuais **não é** a de vários aumentos independentes de longevidade. É a de um **ancestral relativamente longevo com uma redução na linhagem Psittaculidae** (mais ruído de ponta). Isso favorece **H1 (conservação ancestral) + H2 (redução)** e enfraquece, por ora, a premissa de convergência de H3.
- É **um único evento**, o que traz dois limites. Qualquer outra coisa que tenha mudado no mesmo ramo (dispersão para a Australásia/Ásia, dieta nectarívora dos lóris, história reprodutiva) está confundida com ele. E um teste de convergência precisa de ≥ 3 eventos independentes, que não existem aqui (Maddison & FitzJohn 2015).
- A conclusão depende de: árvore (*supertree*), ausência de erro de medida no modelo, um único critério de seleção e cativeiro como contexto. A mesma análise com a árvore de Smith et al. e com erro de medida é obrigatória antes de qualquer afirmação.

## 8. Decisão atualizada

> **Cenário B, reformulado.** O fenótipo contínuo é informativo e agora há **105 espécies com dados genômicos utilizáveis**. A história inferida não sustenta múltiplas origens independentes; sustenta uma **redução em Psittaculidae**.

Consequências para o desenho (registradas no histórico de [`../phase0/analysis_plan.md`](../phase0/analysis_plan.md)):

1. **Teste principal inalterado:** RERconverge/PGLS gene a gene contra o fenótipo contínuo, com massa. Com 105 espécies, a varredura genômica tem poder para efeitos grandes.
2. **Novo teste dirigido:** seleção e relaxamento (aBSREL/RELAX) no **ramo-tronco de Psittaculidae**, que é a hipótese histórica mais bem sustentada. Mas genes com sinal nesse ramo **não** podem ser atribuídos à longevidade sem evidência adicional, porque o ramo é um único evento.
3. **Controle obrigatório:** repetir as associações contínuas **dentro** de Psittaculidae e **dentro** do restante. Uma associação só entre clados é indistinguível do evento único.
4. **Desfecho:** manter a expectativa de vida como principal (fixado antes dos genes), mas **exigir concordância de direção** com a longevidade máxima ajustada para chamar um resultado de robusto. Espécies com `qc_outcomes_inconsistent` entram só na sensibilidade.
5. **Antes de qualquer gene:** auditar os genomas de Hains et al. (BUSCO nos reads/montagens, identidade do espécime por COI/mitogenoma, contaminação). São 96 das 105 espécies.

## 9. Pendências

- Árvore de Smith et al. 2023/2024 (Dryad/Zenodo bloqueados) e repetição da seção 7 com ela.
- Modelo de regimes com erro de medida (p. ex. `l1ou` com `ME`, ou RevBayes).
- Accessions e BUSCO das montagens (NCBI bloqueado). A matriz usa corridas do SRA, não montagens.
- Revisar os recordes do AnAge para as 15 espécies com medidas inconsistentes (fontes primárias/ZIMS).
