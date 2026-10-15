# Sensibilidade à árvore: Smith, Thom & Joseph 2024 (UCE, datada)

15 out. 2026. Sensibilidade obrigatória do plano (`results/phase0/analysis_plan.md` §3). Nenhuma análise de genes foi feita.

## Fonte

Árvore MCC datada (treePL) de Smith, Thom & Joseph 2024, *Bull. Am. Mus. Nat. Hist.* 468: 388 pontas, ~1,14 Mb de UCEs. Arquivo do Zenodo 10867437, espelho do Dryad doi:10.5061/dryad.h70rxwdrx (CC0; o Dryad recusa download por script). md5 `d580615d95bddac463a6bce9f532e1f3`. As árvores bootstrap **datadas** não foram depositadas; o `.ufboot` incluído é a entrada do treePL, sem datação. Por isso só a MCC é usada.

## Harmonização de nomes

`scripts/18_smith_tree.py` → [`smith2024_name_map.tsv`](smith2024_name_map.tsv). Das 385 espécies do grupo interno, 374 entram na árvore final:

- 325 nomes exatos;
- 29 por epíteto único, como *Licmetis* → *Cacatua* e *Palaeornis* → *Psittacula*, todos conferidos;
- 11 manuais, conferidos na Tabela S1. Entre eles, *Clarkona varia* = *Psephotellus varius*, que **não** é um lóri, e dois erros de digitação no arquivo da árvore (*Eunymphicus "cornotus"*, *Myiopsitta "fuchsi"*);
- 11 desdobramentos que a supertree trata como uma espécie só (mantém-se o nome exato);
- 7 mudanças de gênero, 1 sinônimo e 1 grafia.

**As 105 espécies do conjunto genômico estão todas na árvore.**

## Resultados

**Alometria** (PGLS, λ de Pagel), quase idêntica:

| | Supertree (Fase 0b) | Smith 2024 |
| --- | --- | --- |
| Expectativa de vida ~ massa: λ, inclinação | 0,65; 0,35 | 0,66; 0,34 |
| Correlação dos resíduos (longevidade relativa) entre árvores | — | r = 0,9996 (expectativa de vida), 0,988 (máxima) |

**Regimes (PhylogeneticEM, scOU).** Antes, a análise da supertree foi refeita neste computador e reproduziu exatamente a da Fase 0b (K = 3, mesmos deslocamentos).

| | Supertree | Smith 2024 |
| --- | --- | --- |
| K escolhido (LINselect), expectativa de vida | 3: Psittaculidae (101 pontas, −0,35), *Psittaculirostris*, *Tanygnathus* | **1**: só *Psittaculirostris* (−10,8, implausível) |
| Solução com K = 2 e K = 3 | — | Psittaculidae (101 pontas, −0,34) + pontas isoladas; K = 3 **idêntico** à supertree |
| K, longevidade máxima ajustada | 0 | 0 |
| Meia-vida (scOU, K escolhido) | 1,3 Ma | 11,7 Ma |

Assim que o modelo admite um deslocamento num clado, é o de Psittaculidae, com o mesmo tamanho. A seleção de K é instável: na árvore de Smith, o critério prefere um deslocamento numa única ponta extrema. Isso é sintoma de não modelar o erro de medida por espécie, pendência já registrada.

**Teste direto do clado** (PGLS com indicador do clado da Fase 0b; `06_regime_followup.py --clade-from results/phase0b`):

| Desfecho | Supertree | Smith 2024 |
| --- | --- | --- |
| Expectativa de vida relativa | −31 %, p = 0,003 | −25 %, p = 0,052 |
| Longevidade máxima | −23 %, p = 0,038 | −20 %, p = 0,27 |
| Longevidade máxima ajustada | −18 %, p = 0,030 | −17 %, p = 0,12 |
| Sem um subclado de cada vez (expectativa de vida) | −28 % a −35 %, p ≤ 0,015 | −21 % a −28 %, p = 0,018–0,14 |

Detalhes: [`regime_smith2024/shift_clade_contrast.tsv`](regime_smith2024/shift_clade_contrast.tsv).

## Leitura

- **Direção e tamanho** da redução em Psittaculidae são robustos à árvore (−17 % a −31 % conforme desfecho e árvore).
- **A força da evidência não é.** Com a árvore de UCEs, o efeito fica no limiar para a expectativa de vida (p = 0,05) e deixa de ser significativo para a longevidade máxima. O λ maior (0,57 contra 0,40) atribui mais variação à filogenia, e um único evento num clado grande perde poder.
- Pela regra da Fase 0b (mesma direção com a longevidade máxima ajustada), o resultado ainda "concorda em direção", mas não deve ser apresentado como firme. **"Redução em Psittaculidae" passa a ser hipótese de trabalho apoiada, sensível à árvore**, não um achado estabelecido.
- Não refeito com esta árvore: estados ancestrais e pares irmãos. Pendente: modelo de regimes com erro de medida (`l1ou`/RevBayes), que trataria o artefato de ponta isolada.

## Erro de medida por espécie (15 out. 2026)

A expectativa de vida de cada espécie tem erro padrão conhecido (DP posterior de Smeele et al. 2022). O PGLS passou a usar V = σ²·C_λ + diag(SE²) (`fit_pgls_me` em `scripts/lib/phylo.py`, `scripts/19_measurement_error.py`).

**Validação** ([`validation_pgls_me.json`](validation_pgls_me.json)): β, SE e log-verossimilhança idênticos a uma implementação independente em R até < 10⁻⁶; λ e log-verossimilhança idênticos a `phytools::phylosig(method = "lambda", se = )`. O `phytools::pgls.SEy` **não** serve de referência: ele reestima uma escala residual e, na prática, multiplica o SE² conhecido por um fator estimado (outro modelo).

**Teste do clado com erro** ([`measurement_error_contrast.tsv`](measurement_error_contrast.tsv)), expectativa de vida relativa:

| Árvore | Sem erro | Com erro | ΔAIC (com − sem) |
| --- | --- | --- | --- |
| Supertree | −31 %, p = 0,003 | −29 %, p = 0,004 | −18 |
| Smith 2024 | −25 %, p = 0,052 | −23 %, p = 0,073 | −21 |
| Smith, sem Loriini / sem Psittaculini | −28 % / −27 % | −26 % (p = 0,039) / −26 % (p = 0,027) | — |

Modelar o erro melhora muito o ajuste, mas não muda o quadro: o efeito é estável e a evidência depende da árvore.

**As pontas isoladas do PhyloEM são ruído.** As duas espécies que receberam deslocamentos isolados, *Psittaculirostris desmarestii* e *Tanygnathus megalorynchos*, são as duas medidas **mais imprecisas** do conjunto (percentis 99 e 100 do SE). Com o critério fixado antes de olhar (excluir SE > percentil 95; 11 espécies), os regimes foram refeitos:

| Árvore | K escolhido | Melhor deslocamento único (K = 1) | Ganho de log-verossimilhança K0 → K1 |
| --- | --- | --- | --- |
| Supertree | **1: só Psittaculidae** (92 pontas, −0,32) | Psittaculidae | +12,4 |
| Smith 2024 | **0** | Psittaculidae (91 pontas, −0,31) | +6,4 |

O ramo-tronco de Psittaculidae tem o mesmo comprimento nas duas árvores (1,3 Ma; coroa em 36 e 30 Ma), então a diferença vem da topologia e dos ramos dentro e fora do clado.

**Síntese:** em todas as análises, Psittaculidae é o primeiro e o maior deslocamento num clado, sempre com o mesmo tamanho (−0,29 a −0,35 em log; −23 % a −31 %). O suporte estatístico é forte com a supertree e limítrofe com a árvore de UCEs, que é filogeneticamente mais bem fundamentada. Para a H8c e para o desenho dos testes moleculares, Psittaculidae continua sendo o contraste a testar, mas como **hipótese**, sem afirmar que a redução de longevidade está estabelecida.

## Estados ancestrais e pares irmãos (`scripts/21_tree_ancestral_pairs.py`)

| | Supertree | Smith 2024 |
| --- | --- | --- |
| Raiz (BM), longevidade relativa | −0,007 ± 0,34 | −0,010 ± 0,29 |
| Coroa de Psittaculidae | −0,077 ± 0,20 | −0,074 ± 0,17 |
| Pares irmãos entre as espécies com genoma | 29 (8 com Δ ≥ 0,3) | 33 (10 com Δ ≥ 0,3) |
| Pares em comum | 18 | 18 |

Os estados ancestrais não mudam. **Os pares irmãos mudam muito**: 11 pares só existem na supertree e 15 só na árvore de Smith (lista em [`ancestral_pairs_by_tree.json`](ancestral_pairs_by_tree.json); pares da árvore de Smith em [`sister_pairs_smith2024.tsv`](sister_pairs_smith2024.tsv)). Como a árvore de UCEs é a mais bem fundamentada, os contrastes irmãos dos testes moleculares deveriam ser definidos nela, ou restritos aos 18 pares comuns às duas. **Decisão pendente**, a registrar antes de olhar genes (o plano hoje aponta para `phase0b/sister_pair_contrasts.tsv`, da supertree).

## Reproduzir

```bash
.venv/bin/python scripts/18_smith_tree.py
REGIME_TREE=results/phase1_trees/smith2024_species.nwk \
REGIME_MATRIX=results/phase1_trees/species_data_matrix_smith2024.tsv \
REGIME_OUT=results/phase1_trees/regime_smith2024 .renv45/bin/Rscript scripts/R/regime_shifts.R
.venv/bin/python scripts/06_regime_followup.py --tree results/phase1_trees/smith2024_species.nwk \
  --matrix results/phase1_trees/species_data_matrix_smith2024.tsv \
  --out results/phase1_trees/regime_smith2024 --clade-from results/phase0b
```

Erro de medida e pontas ruidosas:

```bash
.venv/bin/python scripts/19_measurement_error.py validate
.venv/bin/python scripts/19_measurement_error.py
.venv/bin/python scripts/19_measurement_error.py se95
REGIME_MATRIX=data/external/matrix_burgio2019_se95.tsv REGIME_OUT=results/phase1_trees/regime_burgio2019_se95 \
  .renv45/bin/Rscript scripts/R/regime_shifts.R
REGIME_TREE=results/phase1_trees/smith2024_species.nwk REGIME_MATRIX=data/external/matrix_smith2024_se95.tsv \
  REGIME_OUT=results/phase1_trees/regime_smith2024_se95 .renv45/bin/Rscript scripts/R/regime_shifts.R
```
