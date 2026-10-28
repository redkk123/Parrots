# Piloto ampliado da H8: 17 genomas de alta qualidade (28 out. 2026)

**Exploratório** (plano, Histórico de 28 out.). São 17 espécies: 10 Psittaculidae (*Melopsittacus*, *Lathamus*, *Psittacula echo* e *P. krameri*, *Agapornis nigrigenis*, *A. fischeri* e *A. roseicollis*, *Loriculus galgulus*, *Neopsephotus*, *Neophema*) × 7 outras (*Strigops*, *Ara ararauna* e *A. rubrogenys*, *Amazona ochrocephala* e *A. vittata*, *Myiopsitta*, *Guaruba*). Não há lóris (nenhuma montagem de qualidade). Foram testados 94 dos 101 genes. Script: `scripts/38_pilot_h8.py` (`PILOT_SET=plus`).

## RELAX (tronco + ramos internos de Psittaculidae × resto; BH no painel, q < 0,10)

| | Genes |
| --- | --- |
| **Relaxados** (k < 1) | **CASP8** (M3, apoptose; k ≈ 0, q = 0,005), **DKC1** (M7; q = 0,004), TXN (M6), NFE2L2 (M5) |
| **Intensificados** (k > 1) | **TERT** (M7; k = 3,2, q = 8 × 10⁻⁶), POT1, ACD (M7), SQSTM1, HSF1, KEAP1 (M5), SOD1, SOD2, CAT (M6) |

- **Heterogeneidade entre módulos:** p = 0,67, **sem sinal de modularidade**.
- O **M7 (telômeros)** é o módulo com mais sinal: 1 relaxado (DKC1) e 3 intensificados (TERT, POT1, ACD).
- M5 e M6 (resposta ao estresse e antioxidantes) têm mais genes intensificados que relaxados.

## aBSREL (seleção episódica; BH entre genes, q < 0,10)

- **Ramo do periquito:** **TERT** (p = 2 × 10⁻⁵), **TERF2** (p = 3 × 10⁻⁴), **RTEL1** (p = 3 × 10⁻¹⁴). **Replica o piloto de 8 espécies.**
- **Tronco de Psittaculidae:** **ACD** (M7; p = 5 × 10⁻⁵).

## Integridade (perda gênica, item 1)

Genes interrompidos com as duas referências (periquito e kakapo) aparecem sobretudo em *Ara rubrogenys* (montagem Nanopore) e *Myiopsitta* (PacBio CLR antigo), as montagens com mais erro de inserção/deleção. O TP53 interrompido aparece em 4 espécies, dentro e fora de Psittaculidae, e o XRCC1 em 7. **Não há perda restrita a Psittaculidae**; esses casos parecem erros de montagem e precisam de reads.

## Leitura (tabela de interpretação pré-registrada)

- **Sem modularidade e sem relaxamento generalizado em Psittaculidae.** Isso enfraquece a versão "perda de travas" da H8 (perda gênica / relaxamento por módulo); H8a/H8b funcionais ficam em aberto.
- **O sinal consistente é a maquinaria de telômeros:** TERT sob seleção intensificada em Psittaculidae e com seleção episódica no ramo do periquito, junto com TERF2 e RTEL1. ACD tem seleção no tronco, DKC1 relaxou. É o candidato mais forte para ligar periquito, câncer e longevidade, e conversa com a literatura de telomerase e câncer.
- **CASP8 relaxado em Psittaculidae** (apoptose, M3) é o único sinal do lado "controle de células danificadas". É compatível com a ideia da H8b (limiar de apoptose), mas é um gene só.
- **Limites:** 17 espécies, sem lóris; montagens não auditadas pelo BUSCO; a análise principal (58 espécies) fica como teste confirmatório.

## Checagem num segundo periquito (Iridian, Illumina; `scripts/39_validate_budgie.py`)

| Gene | VGP | Iridian | Mudanças só do periquito (VGP) | Presentes também no Iridian |
| --- | --- | --- | --- | --- |
| RTEL1 | completo | completo | 63 | **61 (97 %)** |
| POT1 | completo | completo | 26 | **25 (96 %)** |
| MRE11 | interrompido | **completo** | — | — |
| TP53 | interrompido | parcial | — | — |
| TERT, TERF2, TERF2IP | completos | interrompido/parcial | — | — |

- As mudanças de aminoácido exclusivas do periquito nos genes de telômero **se repetem num indivíduo independente, sequenciado com outra tecnologia**. Não são erro da montagem VGP.
- O MRE11 interrompido no VGP está **inteiro** no outro periquito: era erro de montagem.
- TP53 e TERT ficam sem resposta, porque o genoma Iridian é fragmentado demais nesses genes.
