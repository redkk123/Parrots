# Piloto da H8: 101 genes em 8 genomas cromossômicos (28 out. 2026)

**Exploratório** (plano, Histórico de 28 out., 14h50). Com 8 espécies, o resultado não confirma nem refuta a H8; a análise principal (58 espécies) segue rodando. Script: `scripts/38_pilot_h8.py`.

Espécies: 4 Psittaculidae (*Melopsittacus*, *Lathamus*, *Psittacula echo*, *Agapornis nigrigenis*) × 4 outras (*Strigops*, *Ara ararauna*, *Amazona ochrocephala*, *Myiopsitta*). Foram testados 89 dos 101 genes (os demais tinham menos de 4 espécies com ORF completa ou Psittaculidae insuficiente).

## Testes pré-registrados (§10)

- **RELAX** (tronco + ramos internos de Psittaculidae × resto, BH no painel): **nenhum gene relaxado**. Um intensificado: **TERT** (k = 3,7, q = 0,006).
- **Heterogeneidade entre módulos** (média de log2 k, permutação): p = 0,45. Controle (M1–M3) × prevenção/reparo (M4–M6, M8): p = 0,22. **Sem sinal de modularidade** neste piloto.
- **aBSREL, tronco de Psittaculidae:** só RTEL1 (p = 0,017; não sobrevive à correção entre genes).
- **aBSREL, ramo do periquito (H8c):** **TERT** (p = 1,2 × 10⁻⁵), **RTEL1** (p = 1,7 × 10⁻¹⁴), **TERF2** (p = 0,001), **BRCA2** (p = 0,003). Os quatro passam por BH entre os 89 genes a 10 %. Três dos quatro são do **M7 (telômeros)**.

## Descrição complementar (pós-hoc, não pré-registrada)

Fração de aminoácidos em que só o periquito difere das outras 7 espécies ([`budgie_unique_changes.tsv`](budgie_unique_changes.tsv)):

| | Mediana |
| --- | --- |
| M7, telômeros (7 genes) | **5,9 %** |
| Demais módulos (80 genes) | 0,7 % |

Os genes com mais mudanças exclusivas do periquito são RTEL1, TERT e TERF2, seguidos de MCL1, TERF2IP, PMAIP1 e POT1. Cinco dos sete primeiros são de telômero. Nas outras espécies, os mesmos genes ficam em 1–2 %. A diferença entre módulos dá p ≈ 10⁻⁴ por permutação, mas a métrica foi escolhida depois de ver os dados.

## Integridade (item 1)

No genoma de referência do periquito (VGP, GCF_012275295.1), **TP53** e **MRE11** têm 1 base a mais que quebra a leitura. O NCBI marca as duas proteínas como "LOW QUALITY PROTEIN" e as corrige supondo erro de montagem. Pela regra, contam como dado faltante até a checagem com reads.

## Leitura

- Pela tabela de interpretação pré-registrada: sinais restritos a *Melopsittacus* **enfraquecem a H8c** (mudança no tronco de Psittaculidae), sem afetar H8a/H8b.
- **Achado a seguir:** a maquinaria de telômeros (TERT, RTEL1, TERF2, TERF2IP, POT1) mudou muito no ramo do periquito. É compatível com a ligação entre telomerase e câncer, mas aqui não se sabe o efeito funcional. Também pode ser outra coisa: evolução rápida desses genes em aves, ou um ramo longo sem parentes próximos no piloto.
- A análise principal tem lóris (parentes a ~10 Ma do periquito) e diz se o sinal é do periquito ou de um ancestral comum mais recente.
