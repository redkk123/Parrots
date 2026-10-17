# Contexto e Handoff para o Google Colab: Auditoria de Genomas (BUSCO) — Projeto Parrots

Copie e cole as instruções abaixo para o **Gemini no Google Colab**. Ele vai entender todo o projeto e saber exatamente o que precisa ser feito.

---

### Mensagem pronta para colar no Gemini do Colab:

> "Olá, Gemini! Estou rodando uma etapa de bioinformática pesada do meu projeto de pesquisa (genômica comparativa de longevidade em papagaios). 
> 
> **Contexto do Projeto:**
> - Repositório: `https://github.com/redkk123/Parrots.git`
> - Branch de trabalho: `claude/new-session-lwh7xh`
> - Tarefa atual: Rodar a auditoria de qualidade genômica em 96 montagens de genomas de papagaios (Hains et al. 2022 / NCBI).
> - O que o pipeline faz para cada genoma:
>   1. Baixa o FASTA do NCBI Datasets.
>   2. Calcula estatísticas de montagem via `seqkit`.
>   3. Valida a identidade taxonômica pelo genoma mitocondrial (`scripts/11_mito_identity.py`).
>   4. Mede a completude gênica com `compleasm run -l aves --odb odb10`.
>   5. Deleta o FASTA bruto (para economizar espaço) e mantém apenas os relatórios pequenos (`summary.txt`, `full_table.tsv`, `seqkit_stats.tsv`, `time.txt`, `.json`).
> 
> **Estado Atual:**
> - **NÃO precisa começar do zero:** Cerca de 36 dos 96 genomas já foram auditados e commitados no repositório em `results/phase1_audit/busco/` e `results/phase1_audit/mito/`.
> - O script `scripts/10_audit_assembly.sh` **já possui verificação interna**: se o `summary.txt` do genoma já existir, ele pula automaticamente.
> - A fila completa com todos os genomas está em `data/external/audit_queue.tsv` (formato: `GCA_ACCESSION\tSpecies_Name`).
> 
> **Objetivo no Colab:**
> 1. Configurar o ambiente CPU (High-RAM) com `compleasm`, `seqkit`, `ncbi-datasets-cli`, `python 3.10` e dependências.
> 2. Baixar a linhagem `aves_odb10` do compleasm.
> 3. Conectar ao meu Google Drive para ir salvando o progresso de cada genoma em tempo real (garantindo que se o Colab desconectar, nada seja perdido).
> 4. Processar os ~60 genomas restantes em lotes eficientes.
> 
> Por favor, me ajude a estruturar o notebook com as células de execução e monitorar o processo de forma robusta."

---

## Estrutura do Notebook no Colab

Se você quiser rodar direto, o fluxo das células é este:

### Célula 1: Conectar ao Google Drive
```python
from google.colab import drive
import os

drive.mount('/content/drive')
drive_out = '/content/drive/MyDrive/Parrots_results'
os.makedirs(f"{drive_out}/busco", exist_ok=True)
os.makedirs(f"{drive_out}/mito", exist_ok=True)
print("Drive conectado!")
```

### Célula 2: Clonar repositório e instalar ferramentas via Micromamba
```bash
%%bash
cd /content
if [ ! -d "Parrots" ]; then
  git clone -b claude/new-session-lwh7xh https://github.com/redkk123/Parrots.git
fi
cd /content/Parrots
git pull

# Copia a fila de genomas rastreada para a pasta esperada
mkdir -p data/external
cp results/phase1_audit/audit_queue.tsv data/external/audit_queue.tsv

# Baixa micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Cria ambiente com as ferramentas exatas do projeto
./bin/micromamba create -y -p /content/Parrots/.bioenv -c conda-forge -c bioconda \
    compleasm seqkit ncbi-datasets-cli python=3.10 pandas scipy biopython

# Cria atalho .venv para scripts que chamam .venv/bin/python
ln -sfn /content/Parrots/.bioenv /content/Parrots/.venv

# Instala bibliotecas extras do Python
/content/Parrots/.bioenv/bin/pip install pyreadr rdata openpyxl matplotlib

# Baixa a linhagem aves_odb10 requerida pelo protocolo
mkdir -p data/external/busco_lineages
export PATH="/content/Parrots/.bioenv/bin:$PATH"
compleasm download aves --odb odb10 -L data/external/busco_lineages
```

### Célula 3: Execução contínua salvando direto no Drive
```bash
%%bash
cd /content/Parrots
export PATH="/content/Parrots/.bioenv/bin:$PATH"

DRIVE_DIR="/content/drive/MyDrive/Parrots_results"

# Lê a fila e processa apenas os que faltam
while IFS=$'\t' read -r acc sp; do
  # 1. Se já está no repo clonado (os 36 já feitos), pula direto
  if [ -s "results/phase1_audit/busco/$acc/summary.txt" ]; then
    continue
  fi

  # 2. Se já está no Google Drive de execuções anteriores no Colab, pula
  if [ -s "$DRIVE_DIR/busco/$acc/summary.txt" ]; then
    echo "Pulando $acc (já existe no Google Drive)"
    mkdir -p "results/phase1_audit/busco/$acc"
    cp -r "$DRIVE_DIR/busco/$acc"/* "results/phase1_audit/busco/$acc/" 2>/dev/null || true
    continue
  fi

  echo ">>> [$(date '+%H:%M:%S')] Processando $acc ($sp)..."
  bash scripts/10_audit_assembly.sh "$acc" "$sp" 4

  # Copia imediatamente pro Google Drive
  if [ -d "results/phase1_audit/busco/$acc" ]; then
    cp -ru "results/phase1_audit/busco/$acc" "$DRIVE_DIR/busco/"
  fi
  if [ -f "results/phase1_audit/mito/$acc.json" ]; then
    cp -u "results/phase1_audit/mito/$acc.json" "$DRIVE_DIR/mito/"
  fi
  echo ">>> [$(date '+%H:%M:%S')] $acc finalizado e salvo no Drive!"

done < data/external/audit_queue.tsv
```
