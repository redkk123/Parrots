# O que falta fazer

Atualizado em 14 out. 2026. Estado: Fase 0 e Fase 0b concluídas; H8 registrada no protocolo. **Nenhuma análise molecular foi feita ainda.**

Legenda: 🟢 dá para fazer aqui, com o acesso atual · 🔴 precisa liberar rede ou dados · 🧪 precisa de laboratório

---

## 1. Próximo passo (começar por aqui)

### 1.1 🟢/🔴 Auditar os genomas de Hains et al. 2022
São **96 das 105 espécies utilizáveis**. Se o lote tiver problema, o projeto inteiro herda.

- [ ] Montar ou obter montagens e rodar **BUSCO** (aves_odb10). Critério: ≥ 85 % completo no conjunto amplo, ≥ 95 % no painel piloto.
- [ ] **Confirmar a espécie** de cada indivíduo pelo mitogenoma/COI, comparando com a árvore.
- [ ] Checar **contaminação** (humana, bacteriana, entre amostras).
- [ ] Testar se qualidade (BUSCO, N50) se associa à longevidade relativa. Se sim, entra como covariável.
- [ ] Fixar a lista final de espécies antes de olhar qualquer gene.

> 🔴 NCBI, ENA e Dryad estão bloqueados. Sem isso só dá para trabalhar com os reads do SRA (pesado) ou com as montagens do GenomeArk.

---

## 2. Dados e acesso que estão pendentes

- [ ] 🔴 **Liberar rede** (Network access, nas configurações do ambiente): NCBI (Datasets, E-utilities, FTP), Dryad, Zenodo, Crossref, Europe PMC, bioRxiv. Isso destrava quase tudo abaixo.
- [ ] 🔴 **Accessions e versões das montagens** no NCBI (hoje a tabela `assemblies.tsv` tem várias linhas `to_verify`).
- [ ] 🔴 **Árvore de Smith et al. 2023/2024** (UCE, Dryad/Zenodo) e repetir a análise de regimes com ela.
- [ ] 🔴 **Texto completo** das referências (hoje a maioria foi lida só pelo resumo). Prioridade: Ogburn 2001, Wirthlin 2018, Martini 2025, Powers 2019, Kapsetaki 2024.
- [ ] 🔴 **Dados por espécie de Kapsetaki 2024** (necropsias de 108 espécies, com denominadores). É a melhor fonte para o eixo do câncer.
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

- [ ] 🟢 **Modelo de regimes com erro de medida** (`l1ou`, RevBayes). O `PhylogeneticEM` atual não usa o erro padrão de cada espécie.
- [ ] 🟢 **Outras árvores** (assim que existirem): repetir PGLS, estados ancestrais e regimes.
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
- [ ] 🔴 Extrair de Kapsetaki 2024 os valores por espécie e denominadores.
- [ ] 🔴 Comparar *Melopsittacus* com **Loriini** primeiro, depois Psittaculidae (~32 Ma), com *Nymphicus* como controle de manejo.
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

1. **Auditar os genomas de Hains** (seção 1).
2. **Liberar a rede** (seção 2), que destrava NCBI, Dryad e texto completo.
3. **Piloto com 7–9 genomas VGP e ~100 genes** (seção 5.1), que dá para fazer já com o GenomeArk.
4. **Eixo do câncer com dados de necropsia** (seção 7), que é a parte mais forte da H8.
5. **Teste celular H8a × H8b** (seção 6), que exige laboratório e é o único que distingue as duas versões.
