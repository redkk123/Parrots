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

# Checagem com reads: Psittacus erithacus e P. timneh (7 nov. 2026, exploratória)

Montagens: GCA_009867235.2 (*P. erithacus*) e GCA_009867315.2 (*P. timneh*), ambas Iridian com coordenadas idênticas no scaffold (`WYDI02004407.1:1463-9139+` e `WYDJ02004451.1:1463-9139+`). Gene: ACD (M7).

| Parâmetro | Observação | Interpretação |
| --- | --- | --- |
| Coordenada do stop | Posições genômicas 1565–1567 (`TGA`, aa 35 do miniprot, fita forward) no 1º éxon. | Exatamente o mesmo stop identificado em *Loriculus galgulus* (aa 34, 1º éxon). |
| Contexto genômico local | Posições 1500–1650: sequência com **77,3 % de GC** (`GCGCTGGGGGCGCCGTGAGGACGCGGGGGAA...`). | Região hiper-GC / CpG island típica de promotor/éxon 1 aviário. |
| Cobertura por reads curtas | Biblioteca `SRR8236665` (LSU): profundidade média de 31× no corpo do gene ACD (éxons 2–15), mas profundidade **0×** na janela 1500–1650. | Drop-out clássico de sequenciamento Illumina e alinhamento em ilha CpG / extremo 5' divergente. O "stop" do miniprot decorre de alinhamento forçado na borda da montagem desprovida de reads. Não é pseudogenização. |

# Checagem com reads: Myiopsitta monachus (7 nov. 2026, exploratória)

Montagem: GCA_013036005.2 (Iridian, Scaffold). Gene: ACD (M7), coordenada `WUAQ02000083.1:83986751-83994478+`.

| Parâmetro | Observação | Interpretação |
| --- | --- | --- |
| Coordenada do stop | Janela relativa posições 2854–2856 (`TGA`, aa 35 do miniprot, fita forward). | Idêntico a *P. erithacus* (aa 35), *P. timneh* (aa 35) e *L. galgulus* (aa 34). |
| Contexto local | Posições 2750–2900: sequência com **74,7 % de GC** (`AGCTGGAGGCGCCATGAGGACGCGGGGGAGC...`). | Região hiper-GC / CpG de éxon 1 em todas as espécies testadas. |
| Conclusão diagnóstica | O stop do ACD não é perda gênica nem mutação linhagem-específica; é um artefato sistemático de modelagem do miniprot no éxon 1 altamente divergente e rico em GC contra a referência de *Melopsittacus*. | O gene ACD permanece funcional em todas as espécies. Não deve ser considerado interrompido na análise de integridade da H8. |
