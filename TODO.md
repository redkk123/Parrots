# O que falta fazer

Atualizado em 15 out. 2026. Estado: Fase 0 e 0b concluídas; **auditoria das montagens em andamento** (ver [`results/phase1_audit/audit_report.md`](results/phase1_audit/audit_report.md)). **Nenhuma análise de genes foi feita ainda.**

No computador local a rede está liberada (NCBI, ENA, Dryad, Zenodo). Itens 🔴 abaixo passam a ser viáveis.

Legenda: 🟢 dá para fazer aqui, com o acesso atual · 🔴 precisa liberar rede ou dados · 🧪 precisa de laboratório

---

## 1. Próximo passo (começar por aqui)

### 1.1 🟢 Auditar os genomas de Hains et al. 2022 (em andamento)
São **96 das 105 espécies utilizáveis**. Na prática, "Hains" = todas as montagens da Iridian Genomes (235 no NCBI).

- [x] Inventário NCBI de todas as montagens de papagaios (`scripts/09_ncbi_assembly_inventory.py`).
- [ ] ⏳ **BUSCO** (compleasm, aves_odb10) nas 96: lote `scripts/12_audit_batch.sh` rodando (3/96; ~1h20 cada). Critério: ≥ 85 % completo no conjunto amplo, ≥ 95 % no painel piloto.
- [x] **Identidade pelo mitogenoma** (96/96): 61 montagens têm o mitogenoma de **outra** espécie; 22 só têm o de outra ("conflito").
- [ ] ⏳ **Identidade nuclear (Mash)** para resolver os 22 conflitos e medir distância × filogenia (`scripts/16_mash_batch.sh` → `17_mash_identity.py`).
- [ ] **Contaminação nuclear:** medir a fração de reads de outra espécie (ex.: mapear uma amostra de reads contra o mitogenoma próprio e o da doadora; ou k-mers). Saber se a contaminação cruzada vista no mitogenoma chega ao nuclear em nível que importe.
- [x] Qualidade (N50, tamanho, GC) × longevidade relativa: **sem associação** (PGLS). Repetir com BUSCO.
- [ ] **Decidir o critério de N50** (pré-registrado: ≥ 20 kb → só 47 espécies, 12 em Psittaculidae). Manter ou registrar outro **antes** de olhar genes.
- [ ] Rodar BUSCO também nas montagens não-Iridian escolhidas (VGP, B10K etc.) das espécies do conjunto.
- [ ] Fixar a lista final de espécies antes de olhar qualquer gene.

---

## 2. Dados e acesso que estão pendentes

- [x] **Rede liberada** no computador local (NCBI, ENA, Dryad, Zenodo testados em 14 out.).
- [x] **Accessions e versões das montagens**: `assemblies_literature.tsv` resolvida (falta só o lote F1000 de 22 spp.); inventário completo em `results/phase1_audit/assembly_inventory.tsv`.
- [x] **Árvore de Smith et al. 2024** (Zenodo, espelho do Dryad) e regimes repetidos: direção mantida, evidência mais fraca (p = 0,05). Ver `results/phase1_trees/tree_sensitivity_report.md`.
- [x] Estados ancestrais e pares irmãos com a árvore de Smith: ancestrais iguais; **só 18 de 29/33 pares irmãos coincidem**.
- [ ] **Decidir** que pares irmãos usar nos testes moleculares (árvore de Smith ou os 18 comuns) e registrar no plano antes de olhar genes.
- [ ] 🔴 **Texto completo** das referências (hoje a maioria foi lida só pelo resumo). Prioridade: Ogburn 2001, Wirthlin 2018, Martini 2025, Powers 2019, Kapsetaki 2024.
- [x] **Dados por espécie de Kapsetaki 2024** (Zenodo): 15 papagaios com denominadores. Ver `results/phase1_cancer/cancer_report.md`.
- [ ] 🟢 **Revisar os recordes do AnAge** nas 15 espécies com recorde ≈ expectativa de vida (*Ara macao*: 33 anos de recorde contra 35,5 de média). Procurar fonte primária.

---

## 3. Bibliografia e novidade

- [ ] 🔴 **Busca sistemática** (PubMed, Web of Science/Scopus, bioRxiv) com termos registrados em `logs/`. Hoje `novelty_matrix.md` é preliminar e não pode afirmar ineditismo.
- [ ] 🔴 Verificar se **Martini 2025** tem papagaios suficientes para uma análise interna e ver a lista de espécies.
- [ ] 🟢 Preencher os autores marcados `[authors to verify]` em `literature.tsv`.
- [ ] 🔴 Acompanhar o preprint *Parallel selection for longevity* (dez. 2025), se for publicado com escopo maior.
- [ ] 🔴 Procurar séries de **necropsia de lóris** com denominadores (o melhor comparador do periquito).
- [ ] 🔴 Procurar **neuropatologia do envelhecimento** em psitacídeos (nenhuma foi localizada até agora).

---

## 4. Análises ainda da Fase 0

- [x] 🟢 **Erro de medida por espécie** no PGLS (validado contra R) e regimes sem as 5 % pontas mais ruidosas: as pontas isoladas do PhyloEM eram ruído; Psittaculidae continua o primeiro deslocamento nas duas árvores, com suporte forte na supertree e limítrofe na de Smith. Ver `results/phase1_trees/tree_sensitivity_report.md`.
- [ ] Busca de regimes **com** erro por espécie no próprio modelo (RevBayes ou implementação própria do scOU + SE²), se a decisão sobre Psittaculidae precisar de mais que isso.
- [x] 🟢 **Outra árvore** (Smith 2024): PGLS, regimes, ancestrais e pares irmãos refeitos.
- [ ] 🟢 Testar **robustez do achado de Psittaculidae** removendo outros clados (Arini, Cacatuidae) e com outros critérios de seleção.
- [ ] 🟢 **Cobertura estimada → medida:** hoje é bases ÷ 1,2 Gb. Mapear amostras de reads em uma referência para confirmar.

---

## 5. Análise molecular (só depois da seção 1)

Os módulos de genes e os testes estão congelados em `results/phase0/analysis_plan.md`, seção 10. Mudar exige entrada datada.

### 5.1 Extração e alinhamento
- [ ] Extrair os genes dos módulos M1–M8 em cada espécie (miniprot/TOGA com *Melopsittacus* e uma espécie longeva como referência).
- [ ] Ortologia (OrthoFinder ou BUSCO), separar parálogos e multicópia.
- [ ] Alinhamento de códons (MAFFT + PAL2NAL ou MACSE), inspeção de frameshifts e paradas.
- [ ] TERC em análise própria, fora do dN/dS.
- [ ] **Piloto técnico:** 7–9 genomas VGP, ~100 genes, medir tempo e memória (≤ 12 GB RAM) antes de escalar.

### 5.2 Testes principais
- [ ] **RERconverge** contra o fenótipo contínuo (expectativa de vida relativa, ajustada por massa), com permutações filogenéticas e FDR.
- [ ] **Painel dirigido** de ~100 genes como família separada.
- [ ] **RELAX / aBSREL** no tronco de Psittaculidae.
- [ ] Repetir as associações **dentro** de Psittaculidae e **dentro** do restante.
- [ ] Concordância de direção com a **longevidade máxima ajustada** para chamar um resultado de robusto.
- [ ] Leave-one-clade-out (Arini, *Amazona*, Cacatuidae, Loriini, Psittaculini).

### 5.3 Convergência (só se houver suporte)
- [ ] Testar convergência **apenas** se o modelo de regimes sustentar ≥ 3 mudanças independentes. Hoje há um único evento (Psittaculidae).

---

## 6. H8 (travas moleculares)

### Nível 1, genômico (in silico)
- [ ] Perda gênica por **reanotação independente** e cobertura de reads nos módulos M1–M8.
- [ ] Relaxamento da seleção por módulo, com nulo de 1.000 conjuntos aleatórios de genes.
- [ ] Teste de heterogeneidade entre módulos (permutação dos rótulos).
- [ ] Localização: tronco de Psittaculidae × espalhado × só *Melopsittacus*.

### Nível 2, regulatório
- [ ] Procurar dados de **RNA-seq** de papagaios (basal e induzido por estresse) e de ATAC-seq. Se não houver, **registrar a lacuna**.
- [ ] Elementos conservados próximos aos genes dos módulos; variantes de splicing.

### Nível 3, celular 🧪
- [ ] Desenhar o experimento H8a × H8b (dano residual, parada do ciclo, senescência/apoptose, destino das sobreviventes) em fibroblastos de *Melopsittacus*, 1–2 lóris, uma espécie longeva e codorna.
- [ ] Buscar colaboração (laboratório com cultura de células aviárias). Não faz parte do projeto computacional.

---

## 7. Eixo do câncer

- [ ] 🔴 Extrair de Filip & Scope 2002 (texto completo) o **número de animais** e separar benigno de maligno.
- [x] Extrair de Kapsetaki 2024 os valores por espécie e denominadores (rows K24-* em `cancer_comparisons.tsv`).
- [x] Comparar *Melopsittacus* com **Loriini**, Psittaculidae e *Nymphicus*: periquito 20,8 % de neoplasia × Loriini 1,5 %; **mas a calopsita (domesticada) tem 11,4 % e não difere do periquito**. Domesticação é explicação concorrente forte.
- [ ] Separar **domesticação** de linhagem: periquitos não domésticos/selvagens, outra espécie doméstica de Psittaculidae, e **idade à morte** por espécie (pedir aos autores de Kapsetaki 2024).
- [ ] Modelo binomial/beta-binomial hierárquico (espécie, fonte, idade) **só se houver denominadores**.
- [ ] Separar **tecido e malignidade**; considerar *Macrorhabdus*, hormônios, dieta e linhagem doméstica como explicações concorrentes.

## 8. Eixo da neurodegeneração

- [ ] 🔴 Busca específica por cérebros idosos de psitacídeos (zoológicos, serviços de patologia).
- [ ] Preencher `neurodegeneration_evidence_map.tsv` com “não localizado nesta busca” **apenas** quando a busca estiver documentada.
- [ ] Verificar se anticorpos (Aβ, tau, α-sinucleína) foram validados em psitacídeos.

## 9. Eixo oxidativo

- [ ] 🔴 Extrair de Montgomery 2012, Domínguez-de-Barros 2023/2024 e Córdoba-Lanús 2025 os valores, as unidades e a normalização.
- [ ] Não somar ensaios incompatíveis em um índice.
- [ ] Conferir se Montgomery 2012 e os estudos de Tenerife são independentes.

---

## 10. Manuscrito

- [ ] Decidir o título (o plano diz que a convergência só entra se houver suporte).
- [ ] Figuras 1–6 do protocolo. Já existem rascunhos da 1 (árvore), 2 (alometria) e 3 (dois desfechos).
- [ ] Relatar **resultados nulos** e o poder de cada análise.
- [ ] Verificar a novidade antes de escrever a introdução.

---

## 11. Reprodutibilidade (manutenção)

- [ ] Fixar versões dos pacotes Python e R (lockfile).
- [ ] Gravar o resultado do `sha256` dos dados externos a cada execução.
- [ ] Escrever testes para `scripts/lib/phylo.py` (já validado contra R, mas sem testes automáticos).
- [ ] Migrar para Snakemake/Nextflow quando o pipeline molecular entrar.
- [ ] Revisar o que está em `data/processed/` (reads do SRA e metadados).

---

## Resumo por prioridade

1. **Terminar a auditoria** (seção 1.1) e **decidir o critério de N50**.
2. **Separar domesticação no eixo do câncer** (seção 7): o excesso do periquito é grande, mas a calopsita doméstica também é alta. Depois, **texto completo** das referências (seção 2).
3. **Piloto com 7–9 genomas VGP e ~100 genes** (seção 5.1), que dá para fazer já com o GenomeArk.
4. **Eixo do câncer com dados de necropsia** (seção 7), que é a parte mais forte da H8.
5. **Teste celular H8a × H8b** (seção 6), que exige laboratório e é o único que distingue as duas versões.
