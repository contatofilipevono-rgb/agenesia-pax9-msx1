# Prompt para o GPT Codex — finalizar o artigo de agenesia dentária (PAX9/MSX1)

Abra este repositório no Codex. Os arquivos estão na raiz do repositório: onde o texto abaixo diz `agenesia/`, leia a raiz. Trabalhe num branch próprio e entregue as mudanças por pull request. O texto dos artigos de terceiros (PDF e texto do Alibés 2010) não está no repositório; o suplemento do Intarak 2023 está em zip (converta com `scripts/ler_docx_zip.py`, conforme o README).

---

Você vai finalizar um manuscrito científico de odontologia/genética que já tem análise completa, números, figuras e esqueleto escrito. Seu trabalho é **redação final, conferência e empacotamento para submissão**, não refazer a pesquisa.

## Contexto

- **Estudo:** variantes missense de PAX9 (domínio paired) e MSX1 (homeodomínio) na agenesia dentária não sindrômica. Usamos AlphaFold 3 (AlphaFold Server, 41 modelos proteína–DNA), ThermoMPNN (ΔΔG), controle populacional do gnomAD v4 e AlphaGenome 1/2 vs SpliceAI para variantes não codificantes.
- **Achado central:** a maioria das variantes patogênicas do PAX9 age por um de dois mecanismos complementares, contato com o DNA (< 4 Å) ou desestabilização do domínio (ΔΔG ≥ 1,5 kcal/mol). A regra combinada tem sensibilidade 0,84 e especificidade 0,79 fora da amostra (leave-one-out), com LR+ 4,0 e LR− 0,20.
- **Autores:**
  - Filipe Martins Vono (Odontologia, Centro Universitário de Belo Horizonte – UniBH), autor correspondente.
  - Andre Naufel (coautor; afiliação e ORCID a confirmar com o autor — deixe placeholder).
- **Revista-alvo:** European Journal of Oral Sciences (EJOS); plano B: Archives of Oral Biology.
- **Idioma do manuscrito:** inglês.

## Arquivos (pasta `agenesia/`)

- **`PASSAGEM_METODOLOGICA.md`: leia primeiro.** Separa o que foi executado do que foi só planejado. Documenta:
  - a validação do limiar 1,5 e o que ela garante ou não;
  - a independência (ou não) dos conjuntos de calibração, teste e aplicação;
  - as circularidades;
  - os resultados negativos, as limitações e as pendências;
  - as dependências para reproduzir cada número;
  - as correções de 16/11 (Wang 5/8; 16 + 16 patogênicas).

  **Quando divergir de outro arquivo, ela prevalece.**
- `manuscrito_esqueleto.md`: esqueleto em inglês com todas as seções, números e o checklist final (seção 7). **É a base do texto.**
- `LEIAME.md`: diário completo das análises em português, com todos os resultados, decisões e ressalvas. **É a fonte de verdade dos números.**
- `revisao_literatura.md`: revisão de literatura; referências marcadas **[conferir]** ainda não foram verificadas.
- `figuras/`: Fig1–Fig5 (PNG + PDF), Fig1C (estrutura 3D) e `LEGENDAS.md`.
- Tabelas de resultados em CSV:
  - `resultados_afserver.csv`, `grupo_controle.csv`, `validacao_wang2009.csv`, `incertas_classificadas.csv`;
  - `comparacao_ag1_ag2_spliceai.csv`, `splicing_alphagenome_vs_spliceai.csv`, `analise_por_celula.csv`, `expressao_por_celula.csv`;
  - `variantes_PAX9_MSX1.csv`, `literatura_missense.csv`, `jobs_afserver.csv`.
- `scripts/`: código Python que gerou cada número; o uso está no cabeçalho de cada script.
- `estabilidade/`: saídas do ThermoMPNN (posições **0-based** no CSV).
- `estruturas_wt/`: modelos PDB do tipo selvagem.
- `dados_brutos/`: dados das fontes (ClinVar, gnomAD, AlphaMissense, suplemento do Intarak 2023, texto do Alibés 2010).

## Tarefas, nesta ordem

1. **Ler** a passagem metodológica, o esqueleto, o LEIAME e a revisão de literatura inteiros antes de escrever. Os Métodos, a seção de validação e as Limitações devem refletir a passagem: o que está como [PLANEJADO] ou [FALHOU] não pode aparecer como feito.
2. **VUS já reclassificadas** com o limiar validado (contato < 4 Å ou ΔΔG ≥ 1,5): PAX9 45/138 positivas, MSX1 23/72, em `incertas_classificadas.csv`. A seção 3.9 do esqueleto já está atualizada; a versão antiga (ΔΔG ≥ 1) está em `incertas_classificadas_ddg1.csv`. Use os números novos e acrescente na tabela suplementar a coluna de subgrupo (calibração AC ≥ 2 / gnomAD AC = 1 / ausente do gnomAD), como descrito na seção 4.3 da passagem.
3. **Conferir números.** Todo número do texto deve bater com os CSVs ou com o LEIAME. Se algo divergir, recalcule com o script correspondente e registre a correção.
4. **Escrever o manuscrito completo** em inglês, seguindo as normas da EJOS:
   - Consulte as "Author Guidelines" da revista no site da Wiley: limite de palavras, formato do resumo, estilo de referências (provavelmente Vancouver numérico), seções obrigatórias.
   - Siga a estrutura e os parágrafos do esqueleto e expanda a Introdução, os Métodos e a Discussão até os tamanhos indicados.
   - Mantenha o tom cauteloso: resultados "research use only", sem uso clínico.
5. **Referências:**
   - Verifique cada referência marcada [conferir]/[check] em PubMed, Crossref ou Europe PMC (autores, ano, revista, volume, páginas, DOI). Exemplos: Satokata & Maas 1994, Peters 1998, Pejaver 2022, Pak 2023, Nat Commun 2024 sobre homeodomínio e a revisão de Wnt.
   - **Nunca invente referência ou DOI.** Se não conseguir confirmar, deixe `[UNVERIFIED]` e liste no relatório.
   - Inclua as citações das ferramentas: AlphaFold 3 (Abramson 2024), AlphaGenome, ThermoMPNN (Dieckhaus 2024), SpliceAI (Jaganathan 2019), AlphaMissense (Cheng 2023), gnomAD v4, ClinVar, Intarak 2023, Liang 2016, Wang 2009, Alibés 2010.
6. **Pendências do checklist** (seção 7 do esqueleto):
   - **ThermoMPNN:** confirme no artigo/repositório a convenção de sinal (no nosso uso, positivo = desestabilizante) e descreva nos Métodos.
   - **Wang 2009:** confira `validacao_wang2009.csv` contra a Tabela 1 e a Fig. 3 do artigo original, se o acesso for possível. Se não for, marque como pendente.
   - **Ética:** inclua a frase de dados públicos desidentificados, sem participantes humanos.
   - **Termos das ferramentas:** inclua nos Métodos ou nas Declarações que AlphaFold Server e AlphaGenome são de uso não comercial, que saídas do AlphaGenome não têm uso clínico e o aviso exigido pelos AlphaFold Server Output Terms.
7. **Tabelas:** monte as Tables 1–3 e as Suplementares a partir dos CSVs (Table 1 = variantes patogênicas; Table 2 = patogênicas vs gnomAD; Table 3 = concordância com Wang 2009).
8. **Figuras:**
   - Mantenha as figuras e revise as legendas a partir de `figuras/LEGENDAS.md`.
   - Não altere as cores: azul = contato, laranja = desestabilizante, cinza = nenhum/controle, violeta = patogênica. A paleta foi validada para daltonismo.
9. **Declarações:**
   - Data/code availability com placeholder para o DOI do GitHub/Zenodo.
   - Funding: none.
   - Conflitos: none.
   - CRediT: proposta, a confirmar pelos autores.
10. **Cover letter** para a EJOS (1 página).

## Entregáveis (criar pasta `submissao/`)

- `manuscript.docx` (e `manuscript.md`): título, autores com placeholders, resumo estruturado, palavras-chave, texto, referências, legendas.
- `tables.docx`: Tables 1–3.
- `supplementary.xlsx`: VUS reclassificadas, sensibilidade a limiares, AlphaGenome por célula, lista de jobs do AlphaFold Server.
- `cover_letter.docx`.
- As figuras finais renomeadas no padrão da revista (TIFF/PDF na resolução exigida, se as normas pedirem).
- `RELATORIO_FINALIZACAO.md` em português, com:
  - o que mudou;
  - números corrigidos (antes → depois);
  - referências verificadas e as não verificadas;
  - o que ainda depende dos autores (ORCID, afiliação do Andre Naufel, DOI do Zenodo, leitura final).

## Regras

- Não mude conclusões nem números sem recalcular com os scripts; registre toda mudança no relatório.
- Não reescreva a análise do zero nem troque métodos.
- Não há chaves de API no pacote. Os scripts do AlphaGenome precisam de chave e **não devem ser rodados**; use os CSVs já gerados.
- Os modelos do AlphaFold Server (`dados_brutos/afserver/*.zip`) estão no pacote. Os scripts estruturais e estatísticos rodam offline com `python -I -X utf8 scripts/<script>.py ...`, a partir de `agenesia/`; a tabela da seção 10 da passagem lista as entradas de cada um.
- Os scripts do AlphaGenome precisam de chave e **não devem ser rodados**: use os CSVs (há saídas brutas por condição em `*_bruto_*.csv`). Ficou fora só a saída bruta do AlphaGenome 1 (257 MB).
- A integridade metodológica vem antes do impacto. A validação cruzada só mostra que o corte é estável; não chame isso de "validação externa". A reclassificação das VUS é descritiva, porque os controles de calibração estão dentro das VUS. A comparação com Wang 2009 valida o mecanismo, não o classificador.
- Os termos do AlphaFold Server e do AlphaGenome proíbem uso comercial e clínico. O texto não pode sugerir uso diagnóstico.
- Em caso de dúvida científica, não decida sozinho: anote no relatório como pergunta aos autores.
