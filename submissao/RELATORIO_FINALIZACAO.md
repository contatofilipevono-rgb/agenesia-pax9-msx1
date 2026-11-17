# Relatório de Finalização — Manuscrito *PAX9*/*MSX1* (Agenesia Dentária)

**Data:** 16 de novembro de 2026  
**Branch de trabalho:** `feat/manuscript-submission` (em `repo-agenesia/`)  
**Revista-alvo:** *European Journal of Oral Sciences* (EJOS; plano B: *Archives of Oral Biology*)  
**Modelos utilizados:** execução autônoma + duas passagens em `gemini-4-carbon` (`thinking_level="HIGH"`, via `carbon_client.py`: redação científica completa e auditoria adversarial de parecerista), além de verificação programática offline em Python 3.12 (`scipy 1.18.1`) e consultas bibliográficas às APIs do PubMed, Europe PMC e Crossref.

---

## 1. O que foi produzido em `submissao/`

| Arquivo | Conteúdo |
|---|---|
| `submissao/manuscript.md` | Manuscrito completo em inglês: página de título com placeholders, resumo não estruturado (≤ 200 palavras conforme normas da EJOS), 5 palavras-chave MeSH, Introdução (~607 palavras), Material e Métodos (~1.416 palavras), Resultados (~3.771 palavras incluindo Tables 1–3 integradas), Discussão (~1.349 palavras com todas as 12 limitações), Declarações (ética, financiamento, conflitos, CRediT, disponibilidade com avisos obrigatórios do AlphaFold Server e do AlphaGenome), Legendas das Figuras 1–5 (incluindo Fig. 1C) e 27 referências em estilo Vancouver completo. |
| `submissao/manuscript.docx` | Versão Word formatada do manuscrito completo (Times New Roman, espaçamento 1,5, tabelas acadêmicas de três linhas, quebras de página entre folha de rosto, corpo, legendas e referências). |
| `submissao/manuscript.pdf` | Versão PDF limpa para leitura e submissão (folha de rosto, resumo, texto, Tables 1–3, declarações, legendas e referências). |
| `submissao/manuscript_with_figures_preview.pdf` | Versão PDF completa com as Figuras 1–5 (incluindo Fig. 1C) incorporadas ao final para revisão integrada. |
| `submissao/tables.docx` | Arquivo separado com as **Tables 1–3** formatadas em estilo acadêmico de três linhas, notas de rodapé metodológicas e estatísticas consistentes com `manuscript.md`. |
| `submissao/supplementary.xlsx` | Planilha suplementar com 6 abas: `S1_VUS_Classification` (todas as 210 VUS dos domínios com coluna de subgrupo de independência), `S1_VUS_Subgroup_Summary`, `S2_Threshold_Sensitivity` (grade 3×3 completa + ICs de Wilson e Simel + transferência para *MSX1*), `S3_AlphaGenome_By_Cell` (77 variantes não codificantes/splicing por tipo celular + AG1 + SpliceAI), `S3b_Expression_By_Cell` (`expressao_por_celula.csv`) e `S4_AlphaFold_Server_Jobs` (41 jobs + aviso do AlphaFold Server). |
| `submissao/cover_letter.docx` e `cover_letter.pdf` | Carta de submissão de 1 página endereçada ao Editor-in-Chief da *European Journal of Oral Sciences*. |
| `submissao/figuras/` | Cópias em padrão editorial da EJOS: `Figure_1_DNA_distance` a `Figure_5_alphagenome_spliceai` em **TIFF 600 dpi** (compressão LZW) e **PDF vetorial**, mantendo integralmente a paleta validada para daltonismo, além de `Figure_1C_PAX9_AF3_pLDDT.tif` (300 dpi). |
| `submissao/RELATORIO_FINALIZACAO.md` | Este relatório de finalização e auditoria. |

---

## 2. O que mudou e números corrigidos (antes → depois)

### A. Correções herdadas e consolidadas de `PASSAGEM_METODOLOGICA.md`
1. **Concordância com Wang et al. (2009):** **6/8** (corte exploratório antigo ΔΔG ≥ 1,0) → **5/8** no corte principal do artigo (contato < 4 Å ou ΔΔG ≥ 1,5 kcal/mol), pois `S43K` (7,2 Å, ΔΔG = 1,46 kcal/mol) fica 0,04 kcal/mol abaixo de 1,5. O resultado **6/8** com ΔΔG ≥ 1,0 foi mantido explicitamente apenas como análise de sensibilidade.
2. **Composição das patogênicas do *PAX9*:** **19 ClinVar + 13 literatura** (esqueleto §3.1) → **16 ClinVar P/LP + 16 só da literatura (`P-lit`) = 32** (mais *MSX1*: **2 ClinVar + 5 `P-lit` = 7**; total = 39).
3. **VUS dos domínios reclassificadas (ΔΔG ≥ 1,5):** *PAX9* **59/138 → 45/138 (32,6%)**; *MSX1* **24/72 → 23/72 (31,9%)**.
4. **Sobreposição controles × VUS e decomposição em 3 subgrupos:** todos os 124 controles *PAX9* (`AC ≥ 1`) e os 48 controles de calibração (`AC ≥ 2`) estão dentro das 138 VUS; no *MSX1*, 58/59 (`AC ≥ 1`; exceto `S210N`, Likely benign, AC = 143) e 28/29 (`AC ≥ 2`) estão dentro das 72 VUS. Derivada a coluna de subgrupo em `S1_VUS_Classification` e relatadas as taxas separadas:
   - *PAX9*: ausentes do gnomAD **10/14 (71,4%)**; gnomAD `AC = 1` **25/76 (32,9%)**; calibração `AC ≥ 2` **10/48 (20,8%)**.
   - *MSX1*: ausentes do gnomAD **5/14 (35,7%)**; gnomAD `AC = 1` **12/30 (40,0%)**; calibração `AC ≥ 2` **6/28 (21,4%)**.
5. **Enquadramento da validação cruzada:** explicitado em todo o texto (Resumo, Métodos, Resultados, Discussão e Limitações) que a validação *leave-one-out* (80/80 rodadas selecionando 1,5) demonstra **estabilidade do limiar**, e **não** é validação externa; o desempenho fora da amostra coincide com o aparente porque o corte foi idêntico nas 80 dobras; e a comparação com Wang et al. (2009) testa concordância mecanística dentro do próprio conjunto patogênico de calibração.

### B. Novos achados e correções detectados na auditoria deste pacote (scripts + `gemini-4-carbon` HIGH)
1. **Frequência no gnomAD entre as 39 patogênicas (Seção 3.1 vs Table 1):** o rascunho inicial dizia *33 ausentes, 5 singletons (AC = 1) e 1 AC = 2*. A conferência linha a linha na Table 1 e em `variantes_PAX9_MSX1.csv` mostrou que `PAX9 S119F` também tem `AC = 1`, totalizando **32 ausentes (25 *PAX9* + 7 *MSX1*), 6 singletons (`F15I`, `R26W`, `R47W`, `R77G`, `I87F`, `S119F`) e 1 com `AC = 2` (`K91E`)**. Corrigido para **32 / 6 / 1**.
2. **Denominador e percentuais de Alibés et al. (2010):** a auditoria no `gemini-4-carbon` (`HIGH`) apontou inconsistência matemática entre *"88% de 50"* e *"77%"* (`0,77 × 50 = 38,5`). Conferindo diretamente o texto integral de Alibés et al. (2010) (`alibes2010_nar_pax6.txt`), o artigo relata **44 de 50 (88%)** mutações patogênicas explicadas no domínio paired do *PAX6*, sendo **37** envolvendo estabilidade, **6** ligação ao DNA e **9** especificidade de ligação (com sobreposição). Introdução e Discussão foram corrigidas para citar exatamente esses valores.
3. **Valores de p exatos (`scipy.stats`) e ICs de 95% das LRs (pendência §9 da passagem):**
   - Calculados ICs de 95% para as LRs pelo método log-normal de Simel et al. (1991) e por *bootstrap* estratificado (20.000 reamostragens):
     - *PAX9*: `LR+ = 4,05` (Simel `2,29–7,17`; bootstrap `2,50–8,40`); `LR− = 0,20` (`0,197`; Simel `0,09–0,45`; bootstrap `0,04–0,38`); sensibilidade `0,84` (Wilson `0,68–0,93`); especificidade `0,79` (Wilson `0,66–0,88`).
     - *MSX1*: `LR+ = 3,45` (Simel `1,47–8,10`); `LR− = 0,36` (Simel `0,11–1,18`, cruzando 1,0); sensibilidade `0,71` (Wilson `0,36–0,92`); especificidade `0,79` (Wilson `0,62–0,90`).
   - Recalculados todos os testes de Mann–Whitney *U* com `scipy.stats.mannwhitneyu` (com correção de empates e continuidade):
     - Gravidade clínica (*PAX9*): mediana agrupada exata das 25 missense = **10,0** dentes (antes escrita como `≈10`); desestabilizantes vs contato `p = 0,42` (`0,415`, vs `0,40` da aproximação normal sem correção de empates em `gravidade_clinica.py`); nulas vs missense `p = 0,0009` (`0,00088`, arredondado para `0,001` no título da Fig. 4).
     - Testes de Fisher na regra (ΔΔG ≥ 1,5) contra gnomAD `AC ≥ 1` na Table 2 e Seção 3.4: corrigido *PAX9* `1,0×10⁻⁸ → 1,6×10⁻⁸` e *MSX1* `0,041 → 0,045`; e na regra antiga `ΔΔG ≥ 1,0` de *PAX9* (`AC ≥ 2`), `3,0×10⁻⁸ → 3,4×10⁻⁸`.
4. **Valores de Table 3 (Wang et al. 2009):** verificados todos os 8 mutantes diretamente na Tabela 1 e legenda da Figura 3B em `PMC2706687` (HTML completo). Corrigidos os símbolos de EMSA (`G6R: ++++`, `R28P: −/+`, `I87F: −/+`, `K91E: +++`) e as médias de dentes ausentes da Tabela 1 de Wang (`3; 10; 11,3; 13; 6,5; 9; 11,5; 8,4`), que fundamentam a correlação de Spearman `ρ = −0,72` em `validacao_experimental.py`.
5. **Arredondamentos e limites de intervalo:**
   - Mediana de ΔΔG dos 16 variantes *PAX9* de contato: `0,498 kcal/mol` sem arredondamento (`0,50` arredondado; nos scripts com valores já arredondados a 2 casas, a mediana é `0,495`, que em ponto flutuante IEEE-754 arredonda para `0,49` como consta nas legendas antigas). Padronizado para **0,50 kcal/mol** com explicação na Seção 3.3.
   - `pTM` mínimo do AlphaFold 3: `0,84 → 0,83` (`PAX9_S122Y`, `lote2`, em `resultados_afserver.csv`).
   - Expressão em osteoblasto (`expressao_por_celula.csv`): confirmamos que *PAX9* (6,8× *ALB*) e *MSX1* (54,8× *ALB*) também ficam acima do controle negativo em osteoblastos (diferentemente da nota simplificada na passagem §6), e justificamos explicitamente na Seção 2.6 por que `CL:0000134` (célula-tronco mesenquimal) permaneceu como condição primária.
   - Table 2 AlphaMissense: exibido com 4 casas decimais (`0,9999` e `0,9997`) para evitar arredondamento artificial a `1,000`.
6. **Filtragem de variantes não codificantes (84 → 81 → 78 → 77 → 75):** reconciliado no Resumo (`78 → 77` variantes não codificantes/splicing após exclusão da perda do códon de início `g.36662092G>A`) e detalhado integralmente na Seção 2.6 e na legenda da Fig. 5 (`n = 75` SNVs com SpliceAI).
7. **Adequação às normas da *European Journal of Oral Sciences* (EJOS):**
   - Pelas diretrizes da EJOS (Wiley), o resumo deve ser **não estruturado** (parágrafo único, sem subtítulos como *Background/Methods*, máximo **200 palavras**). O resumo em `manuscript.md`/`manuscript.docx` tem **183 palavras por espaço** (**195 tokens alfanuméricos** mesmo se hífens forem divididos pelo portal), garantindo conformidade estrita. *(Nota: caso optem pelo plano B, *Archives of Oral Biology*, será necessário converter o resumo para o formato estruturado.)*
   - Citações numéricas em estilo Vancouver entre parênteses `(1)`, `(1, 2)` em ordem sequencial (`1–27`), com **todos os autores** listados na lista de referências (sem *et al.*) e abreviações NLM/PubMed.
   - Figuras exportadas para `submissao/figuras/` em **TIFF 600 dpi** (e PDF) conforme exigência da EJOS para gráficos/linhas (> 600 dpi).

---

## 3. Verificação de referências e convenção de sinal

### Convenção de sinal do ThermoMPNN (pendência resolvida)
Confirmada diretamente no código-fonte oficial (`Kuhlman-Lab/ThermoMPNN`, `datasets.py`, que inverte os valores `ddG_ML` do conjunto MegaScale: `ddG = -float(row.ddG_ML)`) e no artigo de Dieckhaus et al. (2024, *PNAS* 121(6):e2314853121, `PMC10861915`): **valores positivos de ΔΔG indicam mutações desestabilizantes** ($\Delta\Delta G_{\text{fold}} = \Delta G_{\text{fold,mut}} - \Delta G_{\text{fold,WT}} > 0$), com RMSE de referência de `0,71 kcal/mol` no conjunto MegaScale e `1,51 kcal/mol` no Fireprot.

### Referências verificadas (27/27)
Todas as **27 referências** citadas em `manuscript.md` foram verificadas via APIs do PubMed / Europe PMC / Crossref, incluindo lista completa de autores, volume, fascículo, páginas e DOI:

| # | Referência resumida | PMID / DOI verificado |
|---|---|---|
| 1 | Yu et al. 2019 (*Oral Dis*; revisão via Wnt/agenesia) | PMID `29969831` · `10.1111/odi.12931` |
| 2 | Intarak et al. 2023 (*Jpn Dent Sci Rev*) | PMID `37159578` · `10.1016/j.jdsr.2023.04.001` |
| 3 | Liang et al. 2016 (*Eur J Hum Genet*) | PMID `27381090` · `10.1038/ejhg.2016.78` |
| 4 | Satokata & Maas 1994 (*Nat Genet*) | PMID `7914451` · `10.1038/ng0494-348` |
| 5 | Peters et al. 1998 (*Genes Dev*) | PMID `9732271` · `10.1101/gad.12.17.2735` |
| 6 | Wang et al. 2009 (*Hum Mol Genet*; tabela/legenda conferidas em `PMC2706687`) | PMID `19429910` · `10.1093/hmg/ddp221` |
| 7 | Vastardis et al. 1996 (*Nat Genet*) | PMID `8696335` · `10.1038/ng0896-417` |
| 8 | Yamaguchi et al. 2014 (*PLoS One*) | PMID `25101640` · `10.1371/journal.pone.0102944` |
| 9 | Landrum et al. 2020 (*Nucleic Acids Res*, ClinVar) | PMID `31777943` · `10.1093/nar/gkz972` |
| 10 | Cho et al. 1994 (*Science*, p53–DNA) | PMID `8023157` · `10.1126/science.8023157` |
| 11 | Bullock & Fersht 2001 (*Nat Rev Cancer*, mutantes estruturais/contato de p53) | PMID `11900253` · `10.1038/35094077` |
| 12 | Alibés et al. 2010 (*Nucleic Acids Res*, PAX6) | PMID `20685816` · `10.1093/nar/gkq683` |
| 13 | Schymkowitz et al. 2005 (*Proc Natl Acad Sci U S A*, FoldX) | PMID `16006526` · `10.1073/pnas.0501980102` |
| 14 | Barrera et al. 2016 (*Science*) | PMID `27013732` · `10.1126/science.aad2257` |
| 15 | Kock et al. 2024 (*Nat Commun*, variantes de homeodomínios) | PMID `38600112` · `10.1038/s41467-024-47396-0` |
| 16 | Abramson et al. 2024 (*Nature*, AlphaFold 3) | PMID `38718835` · `10.1038/s41586-024-07487-w` |
| 17 | Dieckhaus et al. 2024 (*Proc Natl Acad Sci U S A*, ThermoMPNN) | PMID `38285937` · `10.1073/pnas.2314853121` |
| 18 | Chen et al. 2024 (*Nature*, gnomAD v4) | PMID `38057664` · `10.1038/s41586-023-06045-0` |
| 19 | Cheng et al. 2023 (*Science*, AlphaMissense) | PMID `37733863` · `10.1126/science.adg7492` |
| 20 | Avsec et al. 2026 (*Nature*, AlphaGenome) | PMID `41606153` · `10.1038/s41586-025-10014-0` |
| 21 | Jaganathan et al. 2019 (*Cell*, SpliceAI) | PMID `30661751` · `10.1016/j.cell.2018.12.015` |
| 22 | Xu et al. 1999 (*Genes Dev*, PDB 6PAX) | PMID `10346815` · `10.1101/gad.13.10.1263` |
| 23 | Hovde et al. 2001 (*Biochemistry*, PDB 1IG7) | PMID `11580277` · `10.1021/bi0108148` |
| 24 | Simel et al. 1991 (*J Clin Epidemiol*, ICs das LRs) | PMID `1941027` · `10.1016/0895-4356(91)90128-v` |
| 25 | Pak et al. 2023 (*PLoS One*) | PMID `36928239` · `10.1371/journal.pone.0282689` |
| 26 | Pejaver et al. 2022 (*Am J Hum Genet*, ClinGen PP3/BP4) | PMID `36413997` · `10.1016/j.ajhg.2022.10.013` |
| 27 | Ramadane-Morchadi et al. 2025 (*Am J Hum Genet*, BRCA1 PP3/BP4) | PMID `40233743` · `10.1016/j.ajhg.2024.12.011` |

*(Nota: nenhuma referência ficou `[UNVERIFIED]`; todas as 27 possuem PMID e DOI conferidos.)*

---

## 4. O que ainda depende exclusivamente dos autores antes do envio

1. **Dados dos coautores e de correspondência:**
   - E-mail institucional e ORCID do autor correspondente **Filipe Martins Vono** (`[CORRESPONDING AUTHOR EMAIL PLACEHOLDER]`, `[ORCID PLACEHOLDER]`).
   - ORCID de **Andre Naufel** (`[ORCID PLACEHOLDER]`, afiliação atualizada para Centro Universitário de Belo Horizonte – UniBH) e confirmação dos papéis CRediT de **Andre Naufel** e **Valeria Gonzaga Santos** (`0009-0009-1333-2518`, American Global Tech University).
2. **Financiamento:** atualizado para registrar apoio parcial por créditos do Google Cloud (sem auxílio monetário externo; acrescentar eventual código de concessão/programa, se desejar especificar).
3. **DOI do Zenodo:** gerar o depósito no Zenodo vinculado ao repositório `https://github.com/contatofilipevono-rgb/agenesia-pax9-msx1` e substituir `[ZENODO DOI PLACEHOLDER]` em `manuscript.md`/`manuscript.docx` e `cover_letter.docx`.
4. **Conferência opcional da subversão do gnomAD (`v4.0` vs `v4.1`) e inspeção visual de Wang 2009:** o endpoint GraphQL `gnomad_r4` não registra a subversão menor no CSV exportado (já declarado como limitação transparente na Seção 4); os valores de Wang 2009 foram conferidos no texto completo HTML em `PMC2706687` (Tabela 1 e legenda da Fig. 3B), restando apenas inspeção visual opcional do gel no PDF original, caso desejem.
5. **Leitura final e revisão de anonimização para a EJOS:** verificar se a revista exigirá envio da página de título em arquivo separado do corpo do manuscrito no portal ScholarOne/Research Exchange da Wiley (para revisão duplo-cega).
