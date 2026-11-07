# Checagem com reads: Loriculus galgulus (7 nov. 2026, exploratória)

Montagem: GCA_046563285.1 (Harvard). Reads: SRR24579577 (HiSeq 2500, Hong Kong; **outro indivíduo**, ~5×), mapeados com minimap2 em janelas de ±3 kb dos lócus. Script: `data/external/readcheck/sra_map.sh` (local).

| Gene | Posição do stop (janela) | Resultado |
| --- | --- | --- |
| PSMD6 | 5897–5899 (TGA, aa 192) | **polimórfico**: ~metade das reads tem C em 5897 (CGA, Arg) e metade T (TGA, stop). Esse indivíduo parece heterozigoto; o stop não é fixo. |
| ACD | 10449–10451 (TGA, aa 34, 1º éxon, identidade 0,50) | **sem informação**: cobertura 1–2 reads no sítio (região de baixa mapeabilidade). Compatível com artefato de alinhamento do 1º éxon, mas não provado. |

Leitura: um stop num alelo, com o outro funcional, não é perda gênica da espécie. Uma espécie, um indivíduo extra, ~5×: serve de sinal, não de prova.

# Checagem com reads: Agapornis personatus (7 nov. 2026, exploratória)

Montagem: GCA_034783155.1 (Iridian, Scaffold). Reads: SRR13483920 (AMNH, 427 Mb) e SRR24579550 (HKU, 2,7 Gb; ~5–7×). Gene: PSMD6 (2 stops prematuros na fita reversa).

| Sítio de stop | Posição (janela) | Evidência em reads | Interpretação |
| --- | --- | --- | --- |
| Stop 1 (aa 293) | 5461–5463 (`TCA` na montagem -> `TGA` fita reversa) | Read SRR13483920 cobre com base `G` na pos 5463 (`CGA` = Arg na fita reversa). Reads HKU têm baixa cobertura exata no trímero (1 read idêntica à montagem). | Polimorfismo / erro de montagem desfeito em indivíduo independente. |
| Stop 2 (aa 342) | 4150–4152 (`TCA` -> `TGA`) | Reads HKU cobrem o éxon adjacente (pos 4189+ com MAPQ 50–60), mas o sítio 4150 tem cobertura muito baixa (2 reads). | Inconclusivo / provável baixa complexidade / artefato local. Não há evidência de fixação de perda em múltiplas linhagens. |
