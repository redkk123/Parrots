# Eixo do câncer: periquito × Loriini em necropsias de zoológico

15 out. 2026. Item 7 do `TODO.md`. Script: `scripts/20_cancer_necropsy.py`.

## Fonte

Kapsetaki et al. 2024, *Evol Med Public Health* 12:105–116, conjunto de dados no Zenodo (10.5281/zenodo.11311870, `Supplementary_data_v2.csv`, CC BY 4.0, md5 `0ebeea809e77949423c90b36e95f78e5`). São 5.729 necropsias de adultos de 108 espécies de aves em 25 instituições da AZA, só espécies com ≥ 20 necropsias. Há **15 espécies de papagaios**, com número de necropsias, prevalência de neoplasia e de malignidade. As contagens foram reconstruídas como prevalência × necropsias (exatas: erro < 10⁻⁷). Os grupos seguem a taxonomia do projeto, porque a coluna de família do conjunto é inconsistente.

## Resultados

| Grupo | Espécies | Necropsias | Neoplasia | Maligna |
| --- | --- | --- | --- | --- |
| **Periquito** (*Melopsittacus*) | 1 | 477 | **20,8 %** (99) | **17,4 %** (83) |
| Loriini (9,75 Ma do periquito) | 6 | 399 | 1,5 % (6) | 0,8 % (3) |
| Outros Psittaculidae (*Agapornis* ×2, *Neophema*) | 3 | 169 | 4,1 % (7) | 1,8 % (3) |
| Fora de Psittaculidae | 5 | 165 | 7,3 % (12) | 4,8 % (8) |
| **Calopsita** (*Nymphicus*, domesticada) | 1 | 70 | **11,4 %** (8) | **10,0 %** (7) |

Testes (periquito contra o grupo; `necropsy_contrasts.tsv`):

| Comparação | Fisher, neoplasia | Fisher, maligna | Beta-binomial entre espécies (log₁₀ p, unilateral) |
| --- | --- | --- | --- |
| × Loriini | OR 17, p ≈ 10⁻²¹ | OR 28, p ≈ 10⁻¹⁹ | −26 / −22 |
| × outros Psittaculidae | OR 6, p ≈ 10⁻⁷ | OR 12, p ≈ 10⁻⁸ | −2,6 / −2,4 |
| × todos os outros papagaios (14 spp.) | OR 7, p ≈ 10⁻²² | OR 11, p ≈ 10⁻²² | **−2,7 / −2,6** |
| × calopsita | OR 2,0, **p = 0,076** | OR 1,9, **p = 0,17** | — |

Entre os outros papagaios (sem o periquito), **domesticados × não domesticados**: neoplasia 8,7 % × 2,5 % (p = 0,004), malignidade 7,8 % × 1,0 % (OR 8,8, p = 0,0002). Os domesticados são a calopsita e *Neophema splendida*.

## Leitura

1. **O excesso do periquito é real e grande** nesta série padronizada (mesmas instituições, necropsia obrigatória, histopatologia). Contra os lóris, seus parentes mais próximos com dados, a diferença é de cerca de 14× em neoplasia e 23× em malignidade. Considerando a variação entre espécies, o periquito fica fora da distribuição dos outros 14 papagaios (p ≈ 0,002).
2. **A domesticação é uma explicação concorrente forte.** A calopsita, de outra família (Cacatuidae, 41 Ma do periquito) e também domesticada, tem malignidade de 10 %, sem diferença significativa do periquito com n = 70. Pela filogenia, a calopsita deveria parecer com os outros papagaios. Ela se parece com o periquito, que compartilha com ela o manejo e a história de criação em cativeiro.
3. **Para a H8:** o contraste periquito × Loriini, o mais forte do eixo, **não separa** "travas moleculares perdidas na linhagem do periquito" de "efeito de linhagens domésticas" (endogamia, seleção para cor e tamanho, obesidade, dieta). Antes de usá-lo como apoio à H8 é preciso de pelo menos um destes:
   - periquitos **selvagens ou de linhagem não doméstica** (populações australianas, coleções de museu);
   - outra espécie doméstica de Psittaculidae (ex.: *Agapornis roseicollis*, *Neophema* em criação intensa) com denominadores;
   - idade à morte (a série tem idade para 1.287 indivíduos, mas os dados públicos não trazem idade por espécie). Periquitos vivem pouco, e câncer cresce com a idade; se os periquitos necropsiados forem relativamente mais velhos (proporção da longevidade), parte do excesso é idade.
4. **Limites:** zoológicos da AZA (não clínica de pets); 15 espécies só; sem tecido de origem nem idade nos dados públicos; *Macrorhabdus* (carcinoma proventricular) não separável.

## Por grupo: papagaios, araras, cacatuas, calopsita

Com o corte de ≥ 20 necropsias, as séries da AZA **não têm** amazonas (*Amazona*), papagaio-cinzento (*Psittacus*) nem araras do gênero *Ara*. O conjunto de 2025 do mesmo grupo (Zenodo 10.5281/zenodo.14834543) lista várias delas, mas só com glicose e longevidade, sem necropsias. Compton et al. 2025 (*Cancer Discovery*) usa o mesmo corte e as mesmas instituições. Por espécie (IC 95 % exato; Fisher contra o periquito; tecido do conjunto de 2025):

| Grupo | Espécie | Necropsias | Neoplasia | IC 95 % | Maligna | p contra periquito | GI / não GI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Periquito | *Melopsittacus undulatus* | 477 | 20,8 % | 17–25 % | 83 | — | 5,9 % / 13,0 % |
| Calopsita | *Nymphicus hollandicus* | 70 | 11,4 % | 5–21 % | 7 | 0,076 | 1,4 % / 8,6 % |
| Cacatua | *Cacatua moluccensis* | 20 | 5,0 % | 0–25 % | 0 | 0,094 | — / 5,0 % |
| Arara | *Anodorhynchus hyacinthinus* (arara-azul) | 31 | 9,7 % | 2–26 % | 1 | 0,17 | 3,2 % / 6,5 % |
| Lóri | *Trichoglossus haematodus* | 190 | 2,6 % | 1–6 % | 3 | 10⁻¹⁰ | 1,6 % / 1,1 % |
| Lóri | *Eos bornea* | 60 | 1,7 % | 0–9 % | 0 | 6 × 10⁻⁵ | 1,7 % / — |
| Agapórnis | *Agapornis fischeri* | 28 | 14,3 % | 4–33 % | 2 | 0,63 | — |
| Agapórnis | *Agapornis nigrigenis* | 108 | 1,9 % | 0–7 % | 0 | 7 × 10⁻⁸ | — |

Leitura:
- Contra cacatua, arara-azul e *A. fischeri*, **não há diferença significativa** espécie a espécie: as amostras são pequenas e os intervalos chegam a 25–33 %. O excesso do periquito só é bem demonstrado contra as espécies com mais necropsias (lóris, *A. nigrigenis*).
- Nos dois domesticados, o excesso está sobretudo **fora do trato gastrointestinal** (13,0 % e 8,6 %). Isso inclui tumores renais, gonadais e de pele, frequentes em periquitos e calopsitas de cativeiro. Os lóris têm pouco dos dois tipos.
- Para amazonas, papagaio-cinzento e *Ara* seria preciso outra fonte: séries de patologia de clínicas de aves ou de zoológicos fora da AZA, com denominadores.

## Outras fontes (busca de 15 out. 2026)

Procuradas séries com **denominador por espécie** para amazonas, papagaio-cinzento, araras e cacatuas. Nenhuma fonte pública traz isso. O que existe:

| Fonte | Tipo | O que dá | Limite |
| --- | --- | --- | --- |
| Garner 2006, *Clinical Avian Medicine* cap. 20 (Northwest ZooPath, 1994–2002) | amostras a um serviço de patologia | Psittaciformes: 220 tumores em 3.545 aves (6,2 %). Tumores por grupo: calopsita 39, amazonas 30, araras 28, cacatuas 25, periquito 25, agapórnis 19, cinzento 8 | **só numeradores** por espécie: dá o perfil de tumores, não a prevalência |
| Megan et al. 2024, *J Avian Med Surg* (banco UC Davis–Reavill, 26.013 amostras, 11.737 necropsias, 1998–2018) | necropsias e biópsias | tumores **adipocíticos** e xantomas na necropsia: 1,3 %. Razão de chances por gênero: *Agapornis* 3,5, *Melopsittacus* 3,4, *Myiopsitta* 2,3, ***Amazona* 1,9**, ***Ara* 0,5** | só tumores de gordura (quase todos benignos, ligados a dieta e idade) |
| Bayley et al. 2026, *Aust Vet J* (papagaio-de-barriga-laranja, *Neophema chrysogaster*, 2010–2019) | necropsias de população de cativeiro para conservação | neoplasia como **causa de morte**: 2 de 641 (0,3 %); adultos 0,9 % | causa de morte, não qualquer tumor; muitos jovens |
| Meza-Madrid et al. 2024, *J Comp Pathol* (UNAM, 252 psitacídeos nativos do México) | necropsias | linfoma multicêntrico como neoplasia mais comum | sem acesso aberto; aves de tráfico, jovens |
| Filip & Scope 2002 (clínica de Viena) | casuística clínica | periquito 20,6 % × outros psitacídeos 3,6 % | outros psitacídeos agrupados; sem acesso aberto |

Para os grupos perguntados:
- **Amazonas:** sem prevalência com denominador. São o 2º grupo com mais tumores no serviço de Garner (30, sobretudo carcinoma espinocelular, pólipos e papilomas cloacais e adenocarcinoma biliar) e têm risco aumentado de tumor de gordura (OR 1,9).
- **Araras:** só a arara-azul tem denominador (3/31, 9,7 %, IC 2–26 %). *Ara* tem risco **menor** de tumor de gordura (OR 0,5). No serviço de Garner, predominam papilomas e pólipos cloacais e adenocarcinoma biliar.
- **Cacatuas:** *C. moluccensis* 1/20 (IC 0–25 %). Em Garner, predominam sarcomas de partes moles.
- **Papagaio-cinzento:** só 8 tumores em Garner; nenhuma prevalência.
- **Calopsita:** 11,4 % (Kapsetaki) e o grupo com mais tumores em Garner (39: carcinoma espinocelular, adenocarcinoma de ovário, sarcomas, seminoma). É o padrão de uma espécie doméstica, como o periquito.

O caso do papagaio-de-barriga-laranja vai na mesma direção da hipótese de domesticação: um Psittaculidae parente de *Neophema splendida*, criado em cativeiro mas **não domesticado**, com câncer raro. A medida é diferente, porém (causa de morte).

**A fonte que resolveria isso é o banco UC Davis–Reavill**: mais de 26 mil amostras de papagaios, com gênero, idade e sexo, já usado para calcular razões de chances por gênero. Pedir aos autores (Beaufrère, UC Davis) a prevalência de neoplasia por gênero, com idade, daria denominadores para amazonas, araras, cacatuas e cinzento de uma vez.

## Próximos passos

- Pedir aos autores (Maley/Boddy) a idade à morte e o tecido por espécie de papagaio, ou buscar a base pública da qual vieram os dados.
- Procurar séries de periquitos não domésticos ou de linhagens de exposição × pet.
- Filip & Scope 2002 (texto completo): separar benigno/maligno e checar a calopsita na série clínica.
