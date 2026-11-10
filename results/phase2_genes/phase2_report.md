# H8, análise principal: 58 espécies (8 nov. 2026)

Conjunto congelado em 28 out. (`results/phase1_audit/analysis_species.tsv`): 58 espécies, 15 delas Psittaculidae. Painel de 101 genes (M1–M8, plano §10), 87 com alinhamento testável no RELAX. Duas árvores: Burgio 2019 (principal) e Smith 2024 (sensibilidade). Números tirados de `h8_relax_by_gene.tsv`, `h8_modules.tsv` e `h8_null_modules.tsv`.

> Revisado em 8 nov.: a primeira versão deste arquivo tinha valores de q errados e conclusões fortes demais. Esta versão substitui aquela.

## 1. RELAX por gene (TEST = tronco + ramos internos de Psittaculidae; BH no painel, q < 0,10)

| Gene | Módulo | k (Burgio / Smith) | q (Burgio / Smith) | Espécies no alinhamento (Psittaculidae) | Observação |
| --- | --- | --- | --- | --- | --- |
| TERT | M7 | 5,5 / 3,1 | 0,0008 / 0,011 | 15 (8) | intensificado; alinhamento com poucas espécies |
| MRE11 | M4 | 2,6 / 2,6 | 0,002 / 0,005 | 49 (10) | intensificado |
| SOD1 | M6 | 32 / 50 | 0,031 / 0,034 | 54 (14) | intensificado; k no limite superior na Smith |
| DKC1 | M7 | 0,49 / 0,49 | 0,031 / 0,060 | 49 (11) | relaxado |
| NQO1 | M5 | 0,34 / 0,33 | 0,024 / 0,023 | 41 (9) | relaxado |
| TSC2 | M2 | ≈0 / ≈0 | 0,024 / 0,022 | 56 (15) | relaxado; k no limite inferior |
| BBC3 | M3 | ≈0 / ≈0 | 3 × 10⁻¹³ / 2 × 10⁻¹³ | 36 (9) | relaxado; k no limite, gene curto (193 códons) e rico em GC |
| RTEL1 | M7 | 0 / 0 | 0,026 / 0,023 | 23 (**4**) | relaxado; só 4 Psittaculidae, frágil |
| CCNE1 | M2 | ≈0 / **2,0** | 0,031 / 0,046 | 49 (12) | **direção oposta nas duas árvores**: não é robusto |

**Atenção com k ≈ 0:** valores de k no limite (7 de 87 genes do painel e 17 de 298 do nulo) costumam indicar pouca informação nos ramos de teste, e não necessariamente relaxamento real. TSC2, BBC3, RTEL1 e CCNE1 estão nesse caso.

## 2. Módulos contra o nulo genômico (298 genes, árvore de Burgio; `37_h8_null.py`)

Estatísticas pré-registradas: proporção de genes relaxados e média de log2 k, com pseudo-painéis pareados por tamanho de alinhamento e taxa de fundo.

| Módulo | n | Prop. relaxada (nulo) | p | Média log2 k (nulo) | p |
| --- | --- | --- | --- | --- | --- |
| M1 | 9 | 0 (0,07) | 1 | −0,10 (−0,71) | 0,59 |
| M2 | 6 | 0,33 (0,06) | 0,09 | −3,69 (−0,49) | **0,024** |
| M3 | 10 | 0,10 (0,04) | 0,66 | −1,37 (0,00) | 0,094 |
| M4 | 11 | 0 (0,09) | 0,81 | −0,05 (−0,79) | 0,38 |
| M5 | 8 | 0,12 (0,06) | 0,80 | +0,35 (−0,68) | 0,19 |
| M6 | 13 | 0 (0,03) | 1 | −0,49 (−0,30) | 0,71 |
| M7 | 8 | 0,25 (0,08) | 0,33 | −1,13 (−0,74) | 0,64 |
| M8 | 22 | 0 (0,03) | 1 | −0,27 (+0,04) | 0,59 |

- O único p < 0,05 é a média de log2 k do M2. **Não sobrevive à correção** para 16 testes (8 módulos × 2 estatísticas; Bonferroni 0,003).
- Ele vem de dois genes com k no limite: TSC2 e CCNE1, e o CCNE1 inverte de direção na árvore de Smith.
- Sensibilidade exploratória, feita depois de ver os resultados (`scratch/null_check.py`): com a **mediana** de log2 k no lugar da média, o M2 cai para p = 0,07 e o M3 fica em p = 0,10. Nenhum módulo chega a p < 0,05.
- Heterogeneidade entre módulos (permutação dos rótulos, proporção relaxada): p = 0,023 na Burgio, mas p = 0,13 na Smith.
- No painel, 23 % dos genes têm RELAX p < 0,05; no nulo, 16 %. A diferença é pequena.

**Leitura:** não há evidência robusta de relaxamento concentrado num módulo. O M2 é no máximo uma tendência, e depende de dois genes com estimativas no limite.

## 3. aBSREL (seleção episódica, BH entre genes, q < 0,10)

- **Ramo do periquito (*Melopsittacus*), nas duas árvores:** RTEL1 (q = 1,5 × 10⁻⁸), MYC (7 × 10⁻⁷), TERF2 (5 × 10⁻⁵), TERT e ACD (0,004–0,005), TERF1 (0,006).
- **Tronco de Psittaculidae:** nenhum gene (menor q = 0,23, RTEL1).
- É o sinal que mais se repete: estava no piloto de 8 espécies, no de 17 e agora no de 58. As mudanças de aminoácido do periquito em RTEL1 e POT1 se confirmaram num segundo indivíduo (`39_validate_budgie.py`).
- **Limite:** o ramo do periquito é longo (espécie de geração curta) e domesticado. Seleção episódica num único ramo terminal não diz se o efeito está ligado à longevidade ou ao câncer.

## 4. Integridade (genes interrompidos)

- 329 chamadas "disrupted" no painel. A maioria é artefato: só aparece com uma das referências (OGG1, MUTYH, MYC, HMOX1, PRDX3 com a proteína do periquito; NQO1 com a do kakapo) ou cai em contigs curtos com a proteína truncada (`disrupted_triage*.tsv`, `disrupted_second_reference.tsv`).
- Seis chamadas sobreviveram aos filtros. Verificação (`readcheck/readcheck_Lgalgulus.md`):
  - **PSMD6, *L. galgulus*:** stop heterozigoto num segundo indivíduo (21 reads C × 20 T). Não é perda fixada.
  - **PSMD6, *A. personatus*:** 1 read de outro indivíduo sem o stop 1; stop 2 sem cobertura. Inconclusivo, tendendo a erro de montagem.
  - **ACD (*Loriculus*, *P. erithacus*, *P. timneh*, *Myiopsitta*):** o stop cai sempre no mesmo códon do 1º éxon, numa região com 75–77 % de GC. Em *P. erithacus*, reads Illumina independentes não cobrem essa região (0×), enquanto o resto do gene tem ~31×. Em *P. timneh* e *Myiopsitta* **não houve checagem com reads**: a conclusão vem só da posição idêntica e do conteúdo GC. Leitura mais provável: artefato de alinhamento do 1º éxon.
- **Conclusão:** nenhum caso de perda gênica confirmada no painel. Isso não prova ausência de perda: TP53 e RB1 não puderam ser avaliados, porque as montagens são fragmentadas nesses genes.

## 5. O que dá para dizer

1. **Perda de travas por módulo (H8, nível 1): sem suporte.** Não há perda gênica confirmada nem relaxamento modular robusto contra o nulo.
2. **Telômeros (M7) são o sinal mais consistente:** TERT intensificado em Psittaculidae nas duas árvores (mas com 15 espécies) e seleção episódica em vários genes de telômero no ramo do periquito, que replica os pilotos. É um candidato a seguir, não uma explicação para o câncer do periquito.
3. **Genes isolados com sinal nas duas árvores e sem k no limite:** MRE11, DKC1, NQO1 e SOD1. São achados de um gene só, dentro de um painel de 87.
4. **Fase 3 concluída (10 nov. 2026):** `RERconverge` nos 7.756 genes do genoma (H6: longevidade relativa, tamanho cerebral relativo e cérebro além da longevidade; H8 intra-clado em Psittaculidae e no restante da ordem; enriquecimento GO BP e os 6 eixos funcionais exploratórios) nas duas árvores com 1.000 permulações filogenéticas. Ver **`results/phase3_genome/phase3_report.md`** e `results/phase3_genome/fig6_phase3_genome_rer.png`.
