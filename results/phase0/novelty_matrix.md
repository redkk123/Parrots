# Matriz de novidade — Fase 0

**Data da busca:** 25 set. 2026. **Fontes consultadas:** busca web geral (resumos, páginas de periódicos, PubMed/PMC via resultados de busca). As bases bibliográficas por API (Crossref, Europe PMC, OpenAlex, NCBI) estavam bloqueadas pela política de rede deste ambiente. Portanto, **esta auditoria não é uma busca sistemática** e deve ser repetida com acesso a essas bases antes de afirmar que o recorte é inédito.

**Nível de verificação:** salvo indicação em contrário, os estudos foram avaliados por resumo e descrições públicas, sem leitura do texto completo (ver `data/curated/literature.tsv`, coluna `verification_level`).

## O que os estudos anteriores fizeram

| Estudo | Escopo taxonômico | Fenótipo | Método de evolução molecular | Trata a história da longevidade *dentro* de Psittaciformes? | Contrastes próximos / incerteza ancestral |
| --- | --- | --- | --- | --- | --- |
| Wirthlin et al. 2018 [L01] | 31 aves; 5 papagaios | longevo vs. demais (categórico) | seleção positiva; elementos regulatórios conservados | Não: papagaios tratados como grupo | Não |
| Martini et al. 2025 [L02] | 141 aves | extremos de longevidade, com controle de massa e filogenia | convergência + rede PPI | Não (toda a classe Aves) | Não específico para papagaios |
| Matsuda & Makino 2024 [L09] | aves vs. morcegos vs. não voadores | voo/metabolismo/longevidade | taxas evolutivas; substituições convergentes | Não | Não |
| *Parallel selection…* 2025, preprint [L10] | mamíferos longevos + aves grandes longevas | longevidade | seleção purificadora e positiva | Não | Não |
| Smeele et al. 2022 [L03] | 244 papagaios (218 com estimativas publicadas) | expectativa de vida contínua (cativeiro) | nenhum (fenotípico) | Sim, mas sem genômica | Modela fenótipo contínuo; não reconstrói histórias |
| Hains et al. 2022 [L07]; 22 spp. F1000 2020 [L08] | 94 + 22 papagaios | — | recurso genômico (short-read) | — | — |

## O que ainda parece aberto (hipótese de lacuna, a confirmar)

1. **Genômica com fenótipo contínuo *dentro* de Psittaciformes.** Nenhum dos estudos localizados associa taxas ou seleção molecular a uma medida contínua de longevidade ajustada por massa ao longo da árvore dos papagaios. Os dados de Smeele et al. (218 espécies) e os ~100+ genomas short-read (L07, L08) tornam esse desenho tecnicamente possível pela primeira vez.
2. **Incerteza ancestral explícita.** Os estudos genômicos anteriores tratam “longevo” como rótulo de ponta. A Fase 0 aqui mostra que a longevidade relativa alta aparece em clados distantes (Arini neotropicais, *Amazona*, *Poicephalus*/*Coracopsis* africanos, *Cacatua*), mas sem mudanças internas bem sustentadas — isso precisa ser modelado, não presumido.
3. **Robustez por linhagem.** Retirar um clado inteiro (por exemplo Arini ou Cacatuidae) e verificar se a associação persiste não foi feito nos estudos acima.
4. **Integração com os três eixos (câncer, neuro, oxidação) no mesmo arcabouço filogenético.** Nenhum estudo localizado cruza o contraste tumoral de *Melopsittacus* com seus parentes filogenéticos reais (Loriini, ~9,75 Ma) e com a genômica de manutenção somática.

## O que seria replicação ou extensão (e deve ser rotulado como tal)

- Reencontrar seleção em TERT, reparo de DNA, controle proliferativo ou antioxidantes em papagaios = **replicação** de L01/L02.
- Usar o painel dirigido de genes candidatos = **teste de hipótese prévia**, não descoberta.
- Correlacionar tamanho cerebral e longevidade = **replicação** de L03.

## Riscos para a novidade

- **Martini et al. 2025** pode já incluir papagaios suficientes para uma análise interna; verificar a lista de espécies no material suplementar.
- **O preprint de 2025 [L10]** pode ser publicado com escopo ampliado; acompanhar.
- **Smith et al. 2024 [L05]** é a taxonomia/árvore mais recente; qualquer revisor esperará seu uso ou justificativa.
- Pode haver trabalhos em andamento com os genomas VGP de papagaios (8 espécies no GenomeArk).

## Pendências para fechar esta matriz

- [ ] Busca sistemática em Web of Science/Scopus/PubMed/bioRxiv com termos fixados (a registrar em `logs/`).
- [ ] Ler os textos completos de L01, L02, L10 e extrair listas de espécies e genes.
- [ ] Verificar autores marcados como “to verify” em `literature.tsv`.
