# Plano de análise (fixado antes da análise molecular)

**Versão:** 1.0, 25 set. 2026. Alterações posteriores devem ser registradas abaixo, em “Histórico”, com data e justificativa, **antes** de rodar os testes que elas afetam.

## 1. Desfecho

- **Principal:** `log_life_expectancy`, a expectativa de vida ao nascer em cativeiro (Smeele et al. 2022), com EP por espécie.
- **Covariável obrigatória:** `log_body_mass`.
- **Forma do fenótipo nas análises moleculares:** o caráter contínuo `log_life_expectancy` com `log_body_mass` como covariável, **ou** a longevidade relativa (resíduo do PGLS λ, inclinação 0,349). Se usar resíduos, propagar a incerteza por reamostragem dos coeficientes do PGLS e do EP das pontas.
- **Sensibilidade 1:** longevidade máxima documentada (coluna separada, a compilar), com o mesmo ajuste por massa.
- **Sensibilidade 2:** excluir espécies com EP(log) > 0,2 (4 espécies).
- **Não usar** “longevo” binário como desfecho principal. Categorias ficam para análises secundárias.

## 2. Covariáveis e justificativa causal

| Variável | Papel | Decisão |
| --- | --- | --- |
| Massa corporal | confundidor (afeta longevidade e taxas moleculares) | incluir sempre |
| Tamanho cerebral relativo | possível mediador ou causa comum | **não** ajustar no modelo principal; análise secundária |
| Tamanho da postura, idade de 1ª reprodução, desenvolvimento | história de vida (mediador/confundidor) | só em sensibilidade, dada a cobertura incompleta |
| Qualidade da montagem (BUSCO, N50 de contig) | confundidor técnico | covariável ou filtro; testar associação com o fenótipo antes |

## 3. Árvores

- **Principal:** *supertree* de Burgio et al. 2019 (a mesma de Smeele et al. 2022).
- **Sensibilidade obrigatória:** árvore datada de Smith et al. 2023/2024 (UCE), após harmonização taxonômica.
- Para taxas gênicas (RERconverge), comprimentos de ramo estimados por gene sobre a topologia de espécies fixada.

## 4. Critérios de inclusão de genomas (fixar antes dos testes)

- BUSCO (aves_odb10) completo ≥ 85 % **e** N50 de contig ≥ 20 kb para o conjunto amplo; ≥ 95 % e escala cromossômica para o painel piloto.
- Um genoma por espécie. Preferir RefSeq ou VGP curado.
- Ausência ou truncamento de gene não conta como perda sem reanotação independente (miniprot/TOGA).
- Testar se BUSCO e N50 se associam à longevidade relativa (PGLS). Se sim, incluí-los como covariáveis.

## 5. Teste principal e família de hipóteses

1. **Principal:** RERconverge, correlação entre taxas relativas por gene e o fenótipo contínuo (com massa), permulações filogenéticas para calibrar p; FDR (BH) em todos os genes elegíveis.
2. **Painel dirigido (família separada):** os ~100 genes da tabela do protocolo, com FDR própria. Declarado como teste de hipótese prévia, não enriquecimento.
3. **Seleção (complementar):** aBSREL/BUSTED nos ramos com longevidade relativa alta **definidos pelo fenótipo** (item 7); RELAX para separar relaxamento de intensificação. Concordância entre ferramentas no mesmo alinhamento não conta como replicação.
4. **Vias:** enriquecimento com background = genes elegíveis testados; controle por tamanho e sobreposição de vias; permutações que preservem a filogenia.

## 6. Robustez obrigatória

- Leave-one-clade-out: Arini, *Amazona*, Cacatuidae, Loriini, Psittaculini.
- Sensibilidade à árvore (item 3) e ao filtro de alinhamento (antes/depois).
- Reportar tamanhos de efeito com IC, não só p.

## 7. Contrastes e ramos de interesse (definidos agora, só pelo fenótipo)

Ramos “longevidade relativa alta” para testes por ramo, **exploratórios**:
- ancestral de Arini + *Pionites*/*Deroptyus* (maior aumento interno estimado: +0,19, EP 0,13);
- ramos terminais das espécies com genoma e longevidade relativa ≥ +0,4: *Ara macao*, *Amazona vittata*, *A. ochrocephala*, *A. aestiva*, *Cacatua galerita*, *C. sanguinea* (se o genoma for confirmado).

Contrastes próximos (pares irmãos ou quase, com genoma nos dois lados, a completar com o lote Hains):
- *Ara macao* (+0,62) × *Ara ararauna* (+0,23);
- *Amazona ochrocephala*/*aestiva* × *Amazona* de longevidade relativa baixa (a identificar entre as espécies com genoma);
- *Melopsittacus* (−0,29) × Loriini (*Psitteuteles goldiei*, −0,16; outros a mapear).

Esses rótulos **não** implicam origens independentes. A convergência só será testada se o modelo de regimes (feasibility §8) sustentar ≥ 3 mudanças independentes.

## 8. Eixos integrados

- **Câncer:** comparar *Melopsittacus* primeiro com Loriini, depois com Psittaculidae (~32 Ma), e usar *Nymphicus* como controle de manejo. Modelo binomial/beta-binomial hierárquico (espécie, fonte) só com denominadores. Separar benigno/maligno e tecido.
- **Neuro:** entregável = matriz de lacunas; nenhuma inferência de proteção.
- **Oxidação:** tabelar por etapa (produção, neutralização, vulnerabilidade, reparo, tolerância). Não combinar ensaios incompatíveis em um índice.

## 9. Recursos

Piloto (7–9 genomas VGP, ~100 genes): CPU, ≤ 12 GB RAM. Medir tempo e memória antes de escalar para ~100 genomas.

## Histórico

- 2026-09-25: v1.0.
