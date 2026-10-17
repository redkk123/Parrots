# Projeto Parrots — regras fixas

- Dono: Angelo (UNICEPLAC, Brasília). Responder em português informal.
- Nenhuma análise de genes antes da lista final de espécies. Módulos M1–M8 e testes estão congelados em
  `results/phase0/analysis_plan.md` §10; qualquer mudança exige entrada datada no Histórico ANTES de olhar genes.
- Não afirmar ineditismo, nem tratar a "redução em Psittaculidae" como estabelecida (depende da árvore).
- Credenciais e códigos de acesso (ex.: Drury-Reavill/UC Davis) NUNCA vão para o repositório, logs versionados
  ou commits. Passar por variável de ambiente (ex.: `DRURY_REAVILL_CODE`).
- Ao usar o banco Drury-Reavill: necropsia = registros com `tissues_min=3` (principal) e `=2` (sensibilidade),
  como orientou H. Beaufrère; citar o banco. O filtro "Post Mortem" não é interpretável.
- E-mails do Angelo: curtos, naturais, sem cara de IA, sem inventar afiliação. Nunca enviar nada sem pedir.
- Máquina: 11 GB RAM / 8 núcleos. Não rodar jobs pesados em paralelo com o lote BUSCO/compleasm (kill 137).
  Rodar coisas leves com `nice`.
- Para matar processos: `pgrep -f '[x]xx'` com colchete; nunca `pkill -f` cru.
- Ambientes: `.venv/` (Python), `.renv45/` (R), `.bioenv/` (micromamba). Commits terminam com `Co-Authored-By`.
