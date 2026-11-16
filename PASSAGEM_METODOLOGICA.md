# Passagem metodológica — agenesia dentária, PAX9/MSX1

*Redigido em 16/11/2026 para quem vai escrever o artigo (GPT Codex) e para os autores (Filipe Vono, Andre Naufel).*

## Como ler este documento

Este documento é a **referência metodológica do artigo**. Se divergir de `manuscrito_esqueleto.md` ou de partes antigas de `LEIAME.md`, **vale este documento**. O LEIAME é um diário cronológico: as seções de 15/11 ainda usam o corte antigo de ΔΔG ≥ 1.

Cada item traz um status:

- **[EXECUTADO]**: rodou, a saída está no pacote e os números foram conferidos em 16/11.
- **[EXECUTADO, COM RESSALVA]**: rodou, mas há uma limitação que precisa ir para o texto.
- **[PLANEJADO]**: foi pensado ou citado no esqueleto, mas **não foi feito**. Não descrever como feito.
- **[FALHOU]**: tentado sem sucesso. Pode aparecer nos Métodos como tentativa, se for útil.

---

## 0. Correções que o artigo precisa absorver (achadas na auditoria de 16/11)

1. **Concordância com Wang 2009: é 5/8 no corte do artigo, e não 6/8.**
   - O 6/8 do esqueleto e do LEIAME foi calculado com ΔΔG ≥ 1.
   - Com a regra do artigo (contato < 4 Å ou ΔΔG ≥ 1,5), a S43K (4–8 Å, ΔΔG 1,46) passa a ser prevista como "mantém a ligação". No experimento ela perde (EMSA "+", reduzida).
   - Erros com 1,5: S43K, G51S e K91E. Com 1,0: G51S e K91E.
   - Relatar 5/8 como principal e citar que a S43K fica no limiar (1,46). Se quiser, mostrar 6/8 com ΔΔG ≥ 1 como sensibilidade.
   - Script: `scripts/validacao_experimental.py` (o corte agora é o 3º argumento, padrão 1,5).
2. **Composição do conjunto patogênico do PAX9: 16 ClinVar (P/LP) + 16 só da literatura (P-lit) = 32.**
   - O esqueleto (seção 3.1) diz 19 + 13. Isso está errado: 19 é o total de P do PAX9 no ClinVar, incluindo fora do domínio.
   - MSX1: 2 ClinVar + 5 P-lit = 7.
   - Fonte: `grupo_controle.csv` × coluna `rotulo` de `variantes_PAX9_MSX1.csv`.
3. **As VUS do domínio e os controles do gnomAD se sobrepõem.** É o problema de independência mais sério; detalhes na seção 4.
   - Os 124 controles do PAX9 (gnomAD AC ≥ 1) são **todos** VUS ou sem classificação no ClinVar (rótulo U), e estão **todos** dentro das 138 VUS reclassificadas.
   - Os 48 controles recorrentes (AC ≥ 2), usados na calibração, também estão nas 138.
   - O mesmo vale para o MSX1: 58 das 59 controles estão nas 72 VUS, e 28 dos 29 recorrentes.
4. **VUS reclassificadas com o corte 1,5** (feito em 16/11): PAX9 45/138 positivas (antes 59); MSX1 23/72 (antes 24).

---

## 1. Dados e curadoria

| Item | Status | Detalhe | Reproduzir |
|---|---|---|---|
| ClinVar PAX9/MSX1 (missense e não codificantes) | [EXECUTADO] | NCBI E-utilities, acesso em 15/11/2026; só registros de um único gene | `scripts/baixar_dados.py`, `scripts/clinvar_nao_codificantes.py` (internet) |
| gnomAD v4 | [EXECUTADO] | API GraphQL, exomas + genomas. **A subversão exata (v4.0/v4.1) não foi registrada**: conferir antes de escrever | `baixar_dados.py` |
| AlphaMissense | [EXECUTADO] | Arquivos por proteína do AlphaFold DB | `baixar_dados.py` |
| Consolidação e rótulos P/B/U | [EXECUTADO] | `variantes_PAX9_MSX1.csv`. P = P/LP no ClinVar sem conflito; B = B/LB ou AF ≥ 0,1%; U = o resto | `scripts/consolidar.py` (offline, com `dados_brutos/`) |
| Literatura (P-lit) | [EXECUTADO] | Intarak 2023 (suplemento, PAX9, 32 missense) e Liang 2016 Tab. 2 (MSX1, 7). Critério P-lit: descrita em família com agenesia, gnomAD AF < 10⁻⁴, sem conflito no ClinVar | `literatura_missense.csv`; fonte em `dados_brutos/literatura/` |
| Exclusões | [EXECUTADO] | Perdas do códon de início (M1V/L/K/I); PAX9 "Asn16Ile" (posição 16 é V na P55771, provável erro de numeração) | — |
| Leitura dos artigos originais de MSX1 G122E/P247S, PAX9 Y143C, A168G | [PLANEJADO] | Estão fora dos domínios, então não afetam as análises principais | — |
| LOVD | [FALHOU] | IP bloqueado (403) em 15/11. **Não foi usado** | — |

Domínios usados: PAX9 paired 6–130 (UniProt 4–130, limitado pela construção 6–137); MSX1 homeodomínio 172–231 (construção 166–237).

---

## 2. Estrutura (AlphaFold 3) e distância ao DNA

- **[EXECUTADO] 41 jobs no AlphaFold Server, 15/11/2026:** 2 tipo selvagem e 39 mutantes.
  - PAX9 6–137 com o DNA de 6PAX (PAX6); MSX1 166–237 com o DNA de 1IG7.
  - Semente automática, 5 amostras por job; usamos a amostra de maior `ranking_score`.
  - Saídas: `dados_brutos/afserver/lote_parcial.zip` (30 jobs) e `lote2.zip` (11); tabela em `resultados_afserver.csv`.
- **[EXECUTADO] Distância ao DNA:** menor distância entre átomos pesados do resíduo e o DNA, no modelo **tipo selvagem**.
  - Classes: contato < 4 Å, próximo 4–8 Å, distante > 8 Å.
  - Linha de base: 33% dos resíduos do domínio paired tocam o DNA.
  - Script: `scripts/distancia_dna.py`.
- **[EXECUTADO] Robustez do contato:** 12 dos 13 sítios patogênicos em contato aparecem em contato nas 5 amostras do modelo selvagem; só a posição 117 varia (`scripts/robustez.py`, parte c).
- **[EXECUTADO, RESULTADO NEGATIVO] A confiança e a interface dos mutantes não discriminam.**
  - ipTM 0,84–0,87 e pTM 0,84–0,89 em todos os modelos, sem colisões.
  - Contatos proteína–DNA nos mutantes do PAX9: −6% a +9%.
  - Mutantes que perdem a ligação no experimento ficam iguais ou com **mais** contatos: L21P +9,1%, R26W +0,3%, R28P +0,6%, I87F +5,8%.
  - Script: `scripts/interface_dna.py`. É um achado metodológico, não uma falha: o AF3 não reproduz a perda de ligação causada por mutação pontual.
- **[EXECUTADO, COM RESSALVA] Enriquecimento de contato contra o fundo do domínio não é significativo:** 13/27 sítios (48% contra 33%), binomial p = 0,074. O enriquecimento só aparece contra o controle populacional (seção 3).
- **[PLANEJADO] AF3 para um grupo de VUS:** não feito.
- **[PLANEJADO] Figura 1B (sítios coloridos por mecanismo na estrutura, PyMOL/ChimeraX):** não feita. Existe só a Fig1C, colorida por pLDDT.

## 3. Estabilidade (ΔΔG) e comparação com a população

- **[EXECUTADO] ThermoMPNN:** github.com/Kuhlman-Lab/ThermoMPNN, commit `370f76e`, pesos `thermoMPNN_default.pt`, torch 2.5.1+cpu.
  - Entrada: cadeia proteica do modelo selvagem do AF3 (`estruturas_wt/`), **sem DNA**.
  - Mutagênese de saturação: 2.508 trocas no PAX9 e 1.368 no MSX1 (`estabilidade/`).
  - Posições **0-based** no CSV; os scripts conferem o aminoácido de referência. Uma primeira rodada com numeração errada foi descartada.
- **[EXECUTADO, COM RESSALVA] Convenção de sinal:** positivo = desestabilizante, **inferido** da mediana positiva de todas as trocas (0,56 e 0,29). **[PLANEJADO]** confirmar no artigo de Dieckhaus et al. 2024 e no README.
- **[EXECUTADO] Mecanismos no PAX9:**
  - ΔΔG mediano: contato 0,49 (n = 16), próximo 1,84 (n = 8), distante 1,81 (n = 8).
  - Não contato contra contato: AUC 0,88, p = 0,0003.
  - Script: `scripts/analise_estabilidade.py`.
- **[EXECUTADO] Patogênicas contra controle gnomAD** (`scripts/grupo_controle.py` → `grupo_controle.csv`).
  - Controle = missense do gnomAD no domínio, sem rótulo P nem P-lit; versões AC ≥ 1 e AC ≥ 2.
  - PAX9 contra recorrentes (n = 48):
    - contato 50% contra 17%, Fisher p = 0,0025;
    - ΔΔG 1,17 contra 0,34, AUC 0,77, p = 4×10⁻⁵;
    - AlphaMissense AUC 0,87.
  - MSX1: ΔΔG p = 0,002 (AUC 0,88); contato sem diferença (29% contra 21%).
  - **Atenção:** o `grupo_controle.py` imprime a regra combinada com ΔΔG ≥ 1 (94% contra 33%). Para o artigo, usar a regra com 1,5: 84% contra 21%, que é a linha "dist<4,0 Å, ΔΔG≥1,5" de `scripts/robustez.py`.
- **[EXECUTADO, COM RESSALVA] O p do Mann-Whitney é por aproximação normal,** implementada à mão, sem correção de empates. Fisher é exato. Se a revista pedir, recalcular com `scipy.stats`.
- **[FALHOU] DDMut** como segundo preditor: "Internal Server Error" até no exemplo da documentação.
- **[PLANEJADO] Segundo preditor de ΔΔG** (DDMut, FoldX, Rosetta): não feito. Hoje o ΔΔG depende de uma única ferramenta.

## 4. Validação do limiar 1,5 e independência dos conjuntos (seção central)

### 4.1 O que foi feito [EXECUTADO]

- **Regra:** a variante é positiva se fica a menos de 4 Å do DNA **ou** tem ΔΔG ≥ c.
- **Dados:** só o PAX9, com 32 patogênicas (P + P-lit) contra 48 variantes recorrentes do gnomAD (AC ≥ 2). São 80 variantes.
- **Validação cruzada deixando uma de fora (leave-one-out):**
  - em cada uma das 80 rodadas, c é escolhido na grade 0; 0,25; …; 3,0 pelo índice de Youden nas outras 79, e a regra é aplicada na variante deixada de fora;
  - o corte escolhido foi **1,5 em todas as 80 rodadas**;
  - desempenho fora da amostra: sensibilidade 0,84, especificidade 0,79, LR+ 4,0, LR− 0,20.
- **Reproduzir:** `python -I -X utf8 scripts/robustez.py grupo_controle.csv dados_brutos/afserver/lote_parcial.zip`. Rodado de novo em 16/11, com saída idêntica.
- **Sensibilidade aos cortes,** sensibilidade/especificidade (mesmo script, parte b):

| Contato | ΔΔG ≥ 0,5 | ΔΔG ≥ 1,0 | ΔΔG ≥ 1,5 |
|---|---|---|---|
| < 3,5 Å | 0,88/0,62 | 0,81/0,71 | 0,72/0,83 |
| < 4,0 Å | 0,97/0,58 | 0,94/0,67 | **0,84/0,79** |
| < 5,0 Å | 0,97/0,54 | 0,97/0,62 | 0,88/0,75 |

### 4.2 O que a validação cruzada garante e o que não garante (escrever isto no artigo)

1. **Só o corte de ΔΔG foi ajustado dentro da validação.** Ficaram fora dela, decididos com acesso aos dados completos:
   - o corte de distância (4 Å, convencional);
   - a forma da regra ("contato OU ΔΔG");
   - a grade de cortes;
   - a escolha do controle (AC ≥ 2).

   O desenho evoluiu ao longo da análise: primeiro contato, depois ΔΔG ≥ 1 (escolhido olhando os dados), depois 1,5. A validação cruzada não corrige essas escolhas de desenho.
2. **Como o corte escolhido foi o mesmo em todas as rodadas,** o desempenho "fora da amostra" é **idêntico** ao desempenho aparente da regra com 1,5 nas 80 variantes. A validação mostra que o corte é **estável**, mas não mede desempenho em dados novos.
3. **Não existe conjunto de teste externo independente para a regra.**
   - As 8 variantes de Wang 2009 estão **todas** no conjunto patogênico usado na calibração.
   - A comparação com Wang testa outro desfecho (perda de ligação no EMSA) nas mesmas variantes. É uma **validação de mecanismo**, não de classificação.
   - Resultado: 5/8 com 1,5; 6/8 com 1,0.
4. **O MSX1 não entrou na validação cruzada,** porque tem só 7 patogênicas. Os valores do MSX1 com 1,5 (sensibilidade 0,71, falso-positivo 0,21, LR+ 3,5, LR− 0,36) são **aparentes, sem validação**, e vêm de um corte escolhido no PAX9. Isso é transferência, não validação.
5. **Os intervalos de confiança não foram calculados.** Com 32 e 48 variantes, os intervalos de 95% são largos Pelo método de Wilson, calculado em 16/11: sensibilidade 27/32 = 0,84 (IC 0,68–0,93); especificidade 38/48 = 0,79 (IC 0,66–0,88). **[PLANEJADO]** Calcular IC das LRs (por exemplo, por bootstrap).

### 4.3 Independência dos conjuntos: o que foi checado em 16/11

| Par de conjuntos | Sobreposição | Consequência |
|---|---|---|
| Patogênicas × controles gnomAD | **0** (disjuntos por construção) | Ok |
| Patogênicas × VUS reclassificadas | **0** | Ok |
| Controles recorrentes (calibração, AC ≥ 2) × VUS reclassificadas | **PAX9 48/48 e MSX1 28/29** estão dentro das VUS | **Não independentes.** As LRs aplicadas às VUS foram estimadas usando parte dessas mesmas VUS como "negativos" |
| Todos os controles (AC ≥ 1) × VUS | PAX9 124/124; MSX1 58/59 | As "138 VUS" são, na prática, 124 variantes populacionais + 14 só do ClinVar |
| Variantes de Wang 2009 × patogênicas de calibração | 8/8 | A validação com Wang não é independente da calibração |

As 138 VUS do PAX9, decompostas (regra com 1,5):

| Subgrupo | n | Positivas |
|---|---|---|
| Só no ClinVar, ausentes do gnomAD | 14 | 10 (71%) |
| gnomAD AC = 1 (controle "todos", fora da calibração) | 76 | 25 (33%) |
| gnomAD AC ≥ 2 (= controle de calibração) | 48 | 10 (21%, igual à taxa de falso-positivo, por construção) |

No MSX1, os mesmos subgrupos são 14 / 30 / 28, com 5 / 12 / 6 positivas.

**Como escrever:**

- O "controle" não é benigno. É um conjunto de variantes populacionais, quase todas VUS ou sem classificação, e pode conter patogênicas raras. Isso puxa o resultado para a hipótese nula, então o teste é conservador.
- A reclassificação das VUS deve ser apresentada como **descritiva**.
- Só as **14 VUS ausentes do gnomAD** e as **76 com AC = 1** ficaram fora da calibração.
- Sugestão: na tabela suplementar, informar o subgrupo de cada VUS, e no texto dar a taxa de positivas nos dois subgrupos fora da calibração. A coluna pode ser derivada cruzando `incertas_classificadas.csv` com `grupo_controle.csv`.

### 4.4 Circularidades conhecidas

- **ClinVar:** algumas classificações P/LP podem ter usado preditores computacionais (critério PP3) e trazem "no assertion criteria" ou um só submissor. Isso afeta os rótulos, não o ThermoMPNN.
- **AlphaMissense contra controle do gnomAD:** o AlphaMissense foi treinado usando frequência populacional como sinal de benignidade. Comparar com variantes do gnomAD é parcialmente circular, e o AUC de 0,87 é otimista. Além disso, ele dá escores altos até aos controles (mediana de 0,94–0,99).
- **P-lit exige gnomAD AF < 10⁻⁴, e o controle exige AC ≥ 1 ou 2:** os grupos diferem em frequência por construção. Isso não afeta distância nem ΔΔG, que não usam frequência, mas afeta o AlphaMissense.
- **ThermoMPNN** foi treinado em dados de estabilidade (megascale), sem rótulos clínicos. Não há circularidade conhecida com os rótulos.
- **Contato com o DNA:** vem de um modelo previsto do PAX9 com o DNA do cristal do PAX6 (6PAX), sem estrutura experimental do PAX9 com DNA.

## 5. Validação externa e fenótipo

- **[EXECUTADO, COM RESSALVA] Wang et al. 2009** (Hum Mol Genet, doi:10.1093/hmg/ddp221).
  - Arquivo `validacao_wang2009.csv`.
  - **Os valores foram extraídos da página do PMC por leitura automática e não foram conferidos** na Tabela 1 e na Fig. 3 do original. O artigo é ambíguo para a R26W: o texto diz "−" e a legenda diz "−/+".
  - "Perde ligação" foi definido como escore de EMSA ≤ 1.
  - Resultado: 5/8 com 1,5 (S43K, G51S e K91E erram) e 6/8 com 1,0. G51S e K91E agem por mecanismos que o modelo não cobre (dominante-negativo e perda de sinergia com MSX1; perda de transativação).
  - Correlações de Spearman com n = 8, só descritivas: ΔΔG × estabilidade modelada por Wang 0,55; ligação × dentes ausentes −0,72.
- **[EXECUTADO] Mecanismo × gravidade** (`scripts/gravidade_clinica.py <intarak2023_supl.txt> resultados_afserver.csv 1.5`).
  - Dentes ausentes, sem terceiros molares, do suplemento de Intarak 2023: um valor por variante.
  - Com 1,5: contato 9,5 (n = 10), desestabilizante 11,0 (n = 11), nenhum 7,0 (n = 4); AUC 0,61, p = 0,40.
  - Spearman ΔΔG × dentes 0,23.
  - Nulas 15,0 (n = 25) contra missense: AUC 0,77, p = 0,001.
  - **Resultado negativo:** o mecanismo não prevê a gravidade.

## 6. Variantes não codificantes (AlphaGenome × SpliceAI)

- **[EXECUTADO] AlphaGenome 2** (pacote `alphagenome` 2.0.0, API v2, janela de 1 Mb, todos os folds).
  - Escores de splicing (sítios, uso, junções) e RNA-seq, restritos ao próprio gene.
  - Exige escolher uma condição celular, e não há tecido dentário.
  - Condição principal: célula-tronco mesenquimal CL:0000134. Sensibilidade: crista neural CL:0011012, língua UBERON:0001723 e osteoblasto UBERON:0002481_CL:0000062.
  - 84 variantes; 78 com escore; 7 P, 37 B e 34 U.
  - Scripts: `scripts/alphagenome_score.py` e `scripts/analise_por_celula.py`.
  - Precisam de chave da API AlphaGenome, que **não está no pacote**; as saídas estão em CSV.
- **[EXECUTADO] AlphaGenome 1** (API v1, todas as trilhas): `scripts/alphagenome1_score.py` e `comparar_versoes.py`. O v1 devolve escores em quantil, que saturam perto de 1; o v2 devolve escore bruto. Não comparar valores absolutos entre as versões.
- **[EXECUTADO] SpliceAI:** escores pré-calculados via Ensembl VEP REST (`scripts/spliceai_vep.py`). O servidor do Broad não respondeu.
- **[EXECUTADO] Escolha da célula:** pela expressão prevista contra ALB (controle negativo) e GAPDH (positivo), em `expressao_por_celula.csv`.
  - PAX9 não é expresso na crista neural (0,7× ALB).
  - A mesenquimal é a única com os dois genes acima do controle negativo (PAX9 11×, MSX1 35×).
  - Na língua, o controle ALB veio alto (1,95), então a razão não serve.
- **Resultados:**
  - Splicing com 6 P contra 37 B: AUC 1,0 no AlphaGenome 2, no AlphaGenome 1 e no SpliceAI (5 P × 36 B).
  - Expressão: AUC 0,94–0,99, conforme a célula.
  - 7 VUS destacadas, robustas em mesenquimal, osteoblasto e língua, e concordantes com o SpliceAI: MSX1 c.470-9G>A e c.470-10C>A; PAX9 c.772-2A>G, c.772-3C>G, c.771+1G>A, c.771+5G>C e c.771+3G>T.
- **Resultados negativos e ressalvas:**
  - O AlphaGenome **não acrescenta** discriminação ao SpliceAI no splicing.
  - Com 6–7 patogênicas, o AUC de 1,0 é frágil, e as variantes canônicas (+1, −2) são casos fáceis.
  - Os achados novos de verdade são os do MSX1 em −9 e −10, que criam um provável novo sítio aceptor.
  - Três variantes da região 3' não traduzida do PAX9 só aparecem na língua: evidência fraca.
  - Não há patogênicas reguladoras confirmadas para validar o diferencial do AlphaGenome em expressão.
- **Termos de uso:** o AlphaGenome é de uso não comercial e suas saídas **não podem ter uso clínico**. Declarar isso nos Métodos.
- **[PLANEJADO] Ensaio de minigene** para as 7 VUS de splicing: só existe como recomendação.
- **[PLANEJADO] Dados de célula de polpa dentária ou germe dentário:** não feito.

## 7. Figuras

- **[EXECUTADO]** Fig1–Fig5 em PNG e PDF, geradas por `scripts/figuras.py`, que já usa o corte de 1,5. As legendas estão em `figuras/LEGENDAS.md`.
- **[EXECUTADO]** Fig1C: captura 3D do complexo PAX9–DNA do AF3, colorida por pLDDT.
- Paleta validada para daltonismo: azul #2a78d6 = contato; laranja #eb6834 = desestabilizante; cinza #8a8984 = nenhum/controle; violeta #6250d6 = patogênica.
- **[PLANEJADO]** Figura com os sítios na estrutura coloridos por mecanismo.

## 8. Limitações (lista para a Discussão)

1. Não há conjunto de teste externo independente para a regra. A validação cruzada só mostra que o corte é estável (seção 4.2).
2. Controles e VUS se sobrepõem; as LRs aplicadas às VUS são otimistas (seção 4.3).
3. O ΔΔG vem de uma única ferramenta (ThermoMPNN, erro típico ≈ 1 kcal/mol) e de uma única estrutura prevista, sem DNA.
4. O contato vem de um modelo previsto, com DNA emprestado do PAX6.
5. O MSX1 tem n = 7, sem validação cruzada.
6. Os valores de Wang 2009 vêm de um único estudo, foram extraídos automaticamente e ainda não foram conferidos.
7. Gravidade: um valor por variante, de famílias com critérios diferentes e com variabilidade entre portadores.
8. Circularidade do ClinVar (PP3) e do AlphaMissense com o gnomAD.
9. Não houve experimento de laboratório.
10. Não há tecido dentário no AlphaGenome.
11. Amostras pequenas, com intervalos largos (sensibilidade 0,68–0,93; especificidade 0,66–0,88).
12. A análise foi exploratória: o desenho mudou depois de ver os dados (corte 1 → 1,5).

## 9. Pendências (antes da submissão)

| Pendência | Quem |
|---|---|
| Conferir as referências marcadas [conferir]/[check] | Codex (PubMed/Crossref), depois os autores |
| Conferir os valores de Wang 2009 no PDF original | Autores ou Codex, se tiver acesso |
| Confirmar a convenção de sinal do ThermoMPNN no artigo de Dieckhaus 2024 | Codex |
| Registrar a subversão do gnomAD (v4.0 ou v4.1) | Autores (refazer a consulta à API se preciso) |
| IC 95% das LRs (os IC de Wilson de sensibilidade e especificidade já estão na seção 4.2) | Codex (offline, com `grupo_controle.csv`) |
| Coluna "subgrupo" nas VUS (calibração / AC = 1 / ausente do gnomAD) | Codex (offline) |
| Figura 1B por mecanismo | Opcional |
| Segundo preditor de ΔΔG | Opcional, reforça |
| Depósito no GitHub + Zenodo (DOI) | Autores |
| Afiliação e ORCID do Andre Naufel; CRediT | Autores |
| Normas da EJOS e frase de ética (dados públicos) | Codex |

## 10. Dependências para reproduzir

**Ambiente usado:**
- Windows 11; Python 3.12.10; numpy 2.5.4; scipy 1.18.1; matplotlib 3.11.2.
- Os scripts de análise usam quase só a biblioteca padrão.
- Rodar com `python -I -X utf8 scripts/<script>.py ...`, a partir da pasta `agenesia/`. O uso de cada script está no cabeçalho.

| Resultado | Script | Entradas | Precisa de |
|---|---|---|---|
| Tabela de variantes | `baixar_dados.py` → `consolidar.py` | APIs | internet |
| Jobs do AF3 | `gerar_jobs.py` → AlphaFold Server (web, cota de 30/dia) | — | conta Google; saídas já em `dados_brutos/afserver/*.zip` |
| `resultados_afserver.csv` | `extrair_afserver.py`, `distancia_dna.py`, `interface_dna.py` | zips do AF Server | offline |
| ΔΔG | ThermoMPNN (externo) | `estruturas_wt/*.pdb` | torch, ThermoMPNN `370f76e`; saídas já em `estabilidade/` |
| Mecanismos | `analise_estabilidade.py` | resultados + estabilidade | offline |
| Controle gnomAD | `grupo_controle.py` | variantes, `lote_parcial.zip`, estabilidade | offline |
| Validação cruzada, sensibilidade, 5 amostras | `robustez.py` | `grupo_controle.csv`, `lote_parcial.zip` | offline |
| VUS | `classificar_incertas.py ... 1.5` | variantes, `grupo_controle.csv`, zip, estabilidade | offline |
| Wang 2009 | `validacao_experimental.py ... 1.5` | `validacao_wang2009.csv`, resultados | offline |
| Gravidade | `gravidade_clinica.py ... 1.5` | `dados_brutos/literatura/intarak2023_supl.txt`, resultados | offline |
| AlphaGenome 1/2 | `alphagenome*_score.py`, `analise_por_celula.py`, `expressao_por_celula.py` | ClinVar não codificantes | **chave AlphaGenome** (não incluída); saídas em CSV |
| SpliceAI | `spliceai_vep.py` | ClinVar não codificantes | internet (Ensembl VEP) |
| Comparação de versões | `comparar_versoes.py` | CSVs do AlphaGenome e SpliceAI | offline |
| Figuras | `figuras.py <pasta agenesia>` | CSVs + estabilidade | matplotlib |

**Arquivos fora do pacote:**
- `alphagenome1_todas_bruto.csv`: 257 MB; enviado à parte como `alphagenome1_todas_bruto.zip`, com 23 MB.
- O `.env` com as chaves.

**Zips do AlphaFold Server reduzidos (16/11).** Para caber no limite de anexo, `dados_brutos/afserver/lote_parcial.zip` e `lote2.zip` contêm só o que os scripts leem, com os mesmos caminhos internos:
- as 5 `summary_confidences` e os 5 modelos `.cif` de cada job;
- o `full_data` da amostra de maior `ranking_score`;
- o `job_request`.

Ficaram de fora: os `full_data` das outras amostras, os templates e os arquivos auxiliares do servidor.

Checagem feita em 16/11, com saída idêntica byte a byte entre os zips completos e os reduzidos:
- scripts rodados: `extrair_afserver.py`, `distancia_dna.py`, `interface_dna.py`, `robustez.py`, `grupo_controle.py` e `classificar_incertas.py ... 1.5`;
- as saídas reproduzem exatamente `resultados_afserver.csv`, `grupo_controle.csv` e `incertas_classificadas.csv` do pacote.

Os zips completos (75 MB) ficam com os autores.
