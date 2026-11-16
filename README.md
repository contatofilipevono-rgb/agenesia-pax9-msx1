# Agenesia dentária: mecanismos estruturais de variantes de PAX9 e MSX1

Análise computacional, apenas para pesquisa, de variantes missense e não codificantes de **PAX9** (domínio paired) e **MSX1** (homeodomínio) na agenesia dentária não sindrômica. Ferramentas usadas: AlphaFold 3 (AlphaFold Server), ThermoMPNN, gnomAD v4, AlphaMissense, AlphaGenome 1/2 e SpliceAI.

**Autores:** Filipe Vono e Andre Naufel — Odontologia, Centro Universitário de Belo Horizonte (UniBH).
**Situação:** manuscrito em preparação (alvo: European Journal of Oral Sciences). Os resultados ainda não foram revisados por pares.

## Por onde começar

| Arquivo | Conteúdo |
|---|---|
| [`PASSAGEM_METODOLOGICA.md`](PASSAGEM_METODOLOGICA.md) | **Referência metodológica.** Separa o que foi executado do que foi só planejado. Traz a validação do limiar, a independência dos conjuntos, os resultados negativos, as limitações, as pendências e o que cada resultado precisa para ser reproduzido. Prevalece sobre os demais arquivos |
| [`PROMPT_CODEX.md`](PROMPT_CODEX.md) | Instruções para a redação final do artigo |
| [`manuscrito_esqueleto.md`](manuscrito_esqueleto.md) | Esqueleto do manuscrito, em inglês |
| [`LEIAME.md`](LEIAME.md) | Diário cronológico das análises (as seções antigas usam o corte ΔΔG ≥ 1; ver a passagem) |
| [`revisao_literatura.md`](revisao_literatura.md) | Revisão de literatura |
| `figuras/` | Figuras 1–5 (PNG/PDF) e legendas |
| `scripts/` | Código de cada análise; o uso está no cabeçalho de cada script |
| `*.csv` | Tabelas de resultados |
| `dados_brutos/` | Dados das fontes públicas e saídas do AlphaFold Server |
| `estabilidade/`, `estruturas_wt/` | ΔΔG do ThermoMPNN e modelos selvagens |

## Reproduzir

Requisitos: Python 3.12; numpy, scipy e matplotlib só para as figuras. Os scripts de análise usam quase só a biblioteca padrão. Rode a partir da raiz do repositório:

```bash
python -I -X utf8 scripts/robustez.py grupo_controle.csv dados_brutos/afserver/lote_parcial.zip
python -I -X utf8 scripts/classificar_incertas.py variantes_PAX9_MSX1.csv grupo_controle.csv dados_brutos/afserver/lote_parcial.zip estabilidade saida.csv 1.5
python -I -X utf8 scripts/validacao_experimental.py validacao_wang2009.csv resultados_afserver.csv 1.5
```

Para a análise de gravidade, gere primeiro o texto do suplemento de Intarak et al. 2023:

```bash
python -I -X utf8 scripts/ler_docx_zip.py dados_brutos/literatura/intarak2023_supl.zip dados_brutos/literatura/intarak2023_supl.txt
python -I -X utf8 scripts/gravidade_clinica.py dados_brutos/literatura/intarak2023_supl.txt resultados_afserver.csv 1.5
```

- Os scripts do AlphaGenome exigem uma chave pessoal da API, que não está no repositório; as saídas já estão nos CSVs.
- A saída bruta do AlphaGenome 1 (257 MB) não está incluída.
- Os zips do AlphaFold Server foram reduzidos ao que os scripts leem; a passagem (seção 10) descreve a checagem de equivalência.

## Termos, avisos e fontes

- **AlphaFold Server:** os modelos em `dados_brutos/afserver/`, os derivados (`estruturas_wt/`, `resultados_afserver.csv`) e a Fig. 1C são saídas do AlphaFold Server, sujeitas aos [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms). Uso apenas não comercial; veja nos termos as demais restrições.
- **AlphaGenome:** saídas obtidas pela API AlphaGenome, de uso não comercial. **Não devem ser usadas para fins clínicos ou diagnósticos.**
- **Uso clínico:** nada neste repositório é ferramenta diagnóstica. A regra estrutural e as razões de verossimilhança são evidência de pesquisa.
- **Fontes de dados públicos:** ClinVar (NCBI), gnomAD v4, AlphaMissense (AlphaFold DB), SpliceAI via Ensembl VEP, PDB 6PAX e 1IG7.
- **Literatura:**
  - suplemento de Intarak et al. 2023, *Jpn Dent Sci Rev* 59:129–137, doi:10.1016/j.jdsr.2023.04.001, redistribuído sem modificação para fins não comerciais, com atribuição;
  - Liang et al. 2016, *Eur J Hum Genet*;
  - Wang et al. 2009, *Hum Mol Genet*, doi:10.1093/hmg/ddp221.

  Textos completos de artigos de terceiros não estão incluídos.
- **Licença do código:** a definir pelos autores. Até lá, todos os direitos reservados.
