# Checagem com reads: Loriculus galgulus (7 nov. 2026, exploratória)

Montagem: GCA_046563285.1 (Harvard). Reads: SRR24579577 (HiSeq 2500, Hong Kong; **outro indivíduo**, ~5×), mapeados com minimap2 em janelas de ±3 kb dos lócus. Script: `data/external/readcheck/sra_map.sh` (local).

| Gene | Posição do stop (janela) | Resultado |
| --- | --- | --- |
| PSMD6 | 5897–5899 (TGA, aa 192) | **polimórfico**: ~metade das reads tem C em 5897 (CGA, Arg) e metade T (TGA, stop). Esse indivíduo parece heterozigoto; o stop não é fixo. |
| ACD | 10449–10451 (TGA, aa 34, 1º éxon, identidade 0,50) | **sem informação**: cobertura 1–2 reads no sítio (região de baixa mapeabilidade). Compatível com artefato de alinhamento do 1º éxon, mas não provado. |

Leitura: um stop num alelo, com o outro funcional, não é perda gênica da espécie. Uma espécie, um indivíduo extra, ~5×: serve de sinal, não de prova.
