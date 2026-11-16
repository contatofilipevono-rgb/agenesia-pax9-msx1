# Revisão de literatura: genética estrutural da agenesia dentária (PAX9 e MSX1)

*Rascunho de 15/11/2026, base para a introdução e a discussão do artigo. As referências marcadas com **[conferir]** foram citadas de memória ou por resumo e precisam ser checadas no original antes de entrar no manuscrito.*

---

## 1. Agenesia dentária não sindrômica

A agenesia dentária é a ausência congênita de um ou mais dentes, sem outra alteração do organismo. É a anomalia mais comum do desenvolvimento dentário. A forma não sindrômica tem forte componente genético.

Em uma revisão de bases de mutações, os autores reuniram cerca de **198 mutações em 15 genes**. Delas, **91,9% estão em sete genes**: AXIN2, EDA, LRP6, MSX1, PAX9, WNT10A e WNT10B (*Genetic analysis: Wnt and other pathways in nonsyndromic tooth agenesis*, PubMed 29969831 **[conferir autores]**).

Os dois fatores de transcrição desse grupo são o **PAX9** e o **MSX1**. Ambos se ligam ao DNA e ligam ou desligam outros genes. São o foco deste estudo.

**PAX9.** A revisão sistemática de **Intarak et al. (2023)** reuniu 40 artigos, 178 pacientes e 61 mutações, sendo 26 que mantêm a proteína inteira (*in-frame*) e 32 nulas.
- O dente mais afetado é o **segundo molar**.
- As trocas de aminoácido no domínio de ligação ao DNA (domínio paired), em especial na região que liga suas duas metades, se associam à ausência do segundo molar.
- As mutações na ponta final da proteína (região C-terminal) se associam mais à ausência do segundo pré-molar.

**MSX1.** **Liang et al. (2016)** revisaram as mutações do MSX1 e mostraram que as trocas de aminoácido associadas à agenesia se concentram no **homeodomínio**, seu domínio de ligação ao DNA. Exemplos: R202P (Vastardis et al., 1996), L211R e T180I (Yamaguchi et al., 2014), A225T, A227E e L230P. Os autores também apontaram que mutações que encurtam a proteína causam agenesia mais grave do que as trocas de um aminoácido.

## 2. Por que PAX9 e MSX1: biologia do desenvolvimento do dente

O dente se forma por uma troca de sinais em sequência entre o epitélio oral e o ectomesênquima, que vem da crista neural da cabeça. As fases são: iniciação, botão, capuz, campânula, formação da coroa e da raiz, e erupção.

- **Onde atuam:** PAX9 e MSX1 são expressos no **mesênquima dentário** nas fases iniciais. Neles, mantêm a produção de sinais como o BMP4, necessários para o dente passar da fase de botão para a de capuz.
- **O que acontece sem eles:** em camundongos sem Pax9 ou sem Msx1, o desenvolvimento do dente para na fase de botão (Satokata & Maas, 1994; Peters et al., 1998 **[conferir]**).
- **Gravidade por dose:** nos humanos, as mutações costumam atingir só uma das duas cópias do gene (heterozigose). O efeito depende da dose. Uma perda parcial gera dentes menores ou de forma alterada. Uma perda maior gera agenesia.

**A sinergia PAX9–MSX1.** As duas proteínas interagem e ativam juntas os promotores de *Bmp4* e *Msx1*. **Wang et al. (2009)** mostraram que a mutação G51S do PAX9 mantém a ligação ao DNA, mas perde essa sinergia com o MSX1 e tem efeito dominante-negativo. Ou seja, nem toda mutação patogênica age por perda de ligação ao DNA.

## 3. Estrutura dos domínios de ligação ao DNA

- **PAX9:** o **domínio paired** (resíduos 4–130) tem duas metades que se ligam ao DNA, ligadas por um trecho de conexão (*linker*). A referência estrutural é a do PAX6, proteína da mesma família, ligado ao DNA (PDB 6PAX; Xu et al., 1999 **[conferir]**).
- **MSX1:** o **homeodomínio** (resíduos 172–231) é formado por três hélices. A terceira, chamada hélice de reconhecimento, se encaixa no sulco maior do DNA. A referência é a estrutura do Msx-1 de camundongo ligado ao DNA (PDB 1IG7; Hovde et al., 2001 **[conferir]**).

## 4. Como mutações de troca de aminoácido inativam fatores de transcrição

A literatura descreve dois mecanismos principais e complementares:

1. **Perda de contato com o DNA (mutações "de contato"):** o resíduo trocado participa diretamente da ligação, e a estrutura da proteína continua intacta.
2. **Desestabilização da proteína (mutações "estruturais"):** a troca desmonta o domínio, e a ligação se perde de forma indireta.

Essa distinção é clássica no **p53**, com seus mutantes de contato e estruturais **[conferir referência]**. Foi testada de forma sistemática na família PAX por **Alibés et al. (2010)**. Os autores usaram o FoldX, um programa de cálculo de energia de proteínas, para modelar as 50 trocas de aminoácido descritas no domínio paired do **PAX6** (causa da aniridia) e chegaram a quatro conclusões:

- **88% das mutações** têm explicação por energia: 77% por estabilidade e 30% por ligação ou especificidade, com sobreposição entre os grupos.
- Os cortes usados foram **0,8 kcal/mol**, que é o erro típico do FoldX, e **1,6 kcal/mol**, o dobro desse erro, para um efeito muito significativo.
- As trocas observadas em espécies próximas (ortólogos), usadas como controle neutro, quase não desestabilizam.
- **Não houve correlação** entre o efeito previsto e a gravidade clínica, com uma exceção: a hipoplasia foveal se associou a mutações que afetam a ligação, e não a estabilidade.

Outros trabalhos ampliam esse quadro:
- **Barrera et al. (2016, *Science*)** testaram variantes de fatores de transcrição em microarranjos de ligação a proteínas. Encontraram 77 variantes em 28 genes que alteram a afinidade ou a especificidade de ligação.
- Um estudo de 2024 na *Nature Communications* testou 92 variantes de homeodomínios **[conferir autores]**. As ferramentas de previsão separaram moderadamente bem as variantes que mudam a afinidade, mas falharam nas que mudam a especificidade.
- **Wang et al. (2009)** são a principal validação experimental no PAX9. Testaram 8 mutantes do domínio paired (G6R, L21P, R26W, R28P, S43K, G51S, I87F, K91E) por EMSA, um ensaio de ligação ao DNA, e por ensaios de ativação de promotores.
  - Em 7 dos 8, a perda de ligação e de ativação acompanhou a gravidade clínica.
  - A exceção foi a G51S, pelo mecanismo descrito na seção 2.

**Lacuna:** não existe uma análise sistemática do tipo de Alibés para **PAX9 e MSX1 na agenesia dentária**. Também não existe nenhuma análise que use as ferramentas atuais: AlphaFold 3, preditores de estabilidade treinados em grandes bases de dados e controles de população do gnomAD.

## 5. Ferramentas computacionais atuais

| Ferramenta | O que faz | Limitação relevante |
|---|---|---|
| **AlphaFold 3** (Abramson et al., 2024) | Prevê a estrutura de complexos de proteína com DNA, RNA e ligantes | Pouco sensível a mutações pontuais. Tende a manter o mesmo dobramento e a mesma interface (Pak et al., 2023 **[conferir]**). Nossos dados confirmam isso: o ipTM e os contatos com o DNA não mudam com as mutações |
| **ThermoMPNN** (Dieckhaus et al., 2024 **[conferir]**) | Estima o ΔΔG de enovelamento a partir da estrutura, com treino no conjunto *mega-scale* | Erro típico em torno de 1 kcal/mol. Não estima a ligação ao DNA |
| **FoldX** e **DDMut** | Estimam o ΔΔG por campo de força (FoldX) ou por aprendizado de máquina (DDMut) | O FoldX exige licença acadêmica. O DDMut estava fora do ar em 15/11/2026 |
| **AlphaMissense** (Cheng et al., 2023 **[conferir]**) | Dá um escore de patogenicidade para cada troca de aminoácido do genoma | Satura em domínios muito conservados, onde quase todas as trocas recebem escore alto |
| **SpliceAI** (Jaganathan et al., 2019 **[conferir]**) | Prevê alterações de splicing | Só avalia splicing |
| **AlphaGenome 1 e 2** (Avsec et al., 2026) | Prevê o efeito de variantes na expressão, no splicing, na cromatina etc., em janelas de até 1 Mb. A versão 2 é específica por tipo de célula | Não há tecido dentário entre as cerca de 2.800 condições. Muitas trilhas são estimadas pelo modelo, não medidas |

**Integração clínica.** As diretrizes do ACMG/AMP usam evidência computacional nos critérios PP3 (a favor de patogenicidade) e BP4 (a favor de benignidade). A ClinGen recomenda **calibrar** essa evidência em razões de verossimilhança (Pejaver et al., 2022 **[conferir]**). Estudos recentes mostram que informação estrutural (ΔΔG, acessibilidade ao solvente) melhora a calibração. Um exemplo é o BRCA1, em que a estrutura aumentou a precisão da evidência PP3/BP4 (estudo em PMC12120176).

## 6. O que este estudo acrescenta (posicionamento)

1. **Primeira análise sistemática de estrutura e estabilidade para PAX9 e MSX1 na agenesia.** São 39 trocas de aminoácido patogênicas, do ClinVar e da literatura. Cada uma foi classificada por contato com o DNA no AlphaFold 3 e por ΔΔG no ThermoMPNN.
2. **Controle de população**, mais forte que o controle por ortólogos de Alibés. Usamos as variantes do gnomAD v4 dentro do domínio.
   - No PAX9, com a regra do artigo (contato ou ΔΔG ≥ 1,5), contato ou desestabilização aparecem em 84% das patogênicas, contra 21% do controle recorrente. Com ΔΔG ≥ 1, os valores são 94% contra 33% (p = 3×10⁻⁸).
   - Na validação cruzada, a regra tem sensibilidade de 0,84, especificidade de 0,79, LR+ de 4,0 e LR− de 0,20.
3. **O corte de ΔΔG convergiu de forma independente:** a validação cruzada escolheu **1,5 kcal/mol**, perto do corte de 1,6 de Alibés et al.
4. **Validação experimental:** com ΔΔG ≥ 1,5, a previsão de perda de ligação bate em 5 de 8 mutantes de Wang et al. (2009); com ΔΔG ≥ 1, em 6 de 8.
   - Erros: S43K (ΔΔG 1,46, logo abaixo do corte), G51S e K91E. As duas últimas agem por outros mecanismos já descritos.
   - As 8 variantes fazem parte do conjunto de calibração.
5. **Achado metodológico:** o AlphaFold 3 **não reproduz** a perda de ligação causada por uma única mutação. O efeito precisa vir da posição do resíduo e do ΔΔG.
6. **Gravidade:** assim como no PAX6, o mecanismo previsto **não prevê** o número de dentes ausentes. As trocas de aminoácido são menos graves que as mutações nulas (cerca de 10 contra 15 dentes, p = 0,001), o que sugere função residual.
7. **Variantes fora da região codificante:** pelo que pesquisamos, é a primeira aplicação do **AlphaGenome** em odontologia.
   - Ele concordou com o SpliceAI (AUC de 1,0 nas duas ferramentas) e apontou 7 variantes incertas de splicing.
   - As duas mais fortes ficam no íntron do MSX1: c.470-9G>A e c.470-10C>A.
8. **Produto prático:** uma classificação das 210 variantes incertas dos dois domínios, com razões de verossimilhança, para pesquisa e futura curadoria.

## 7. Lições do artigo de referência para o nosso manuscrito

Do artigo de Alibés et al. (2010):
- **Justificar os cortes pelo erro da ferramenta:** 0,8 e 1,6 kcal/mol para o FoldX. Para o ThermoMPNN, citar o erro relatado pelos autores e mostrar que o corte de 1,5 escolhido pela validação cruzada é coerente com esse erro.
- **Manter a estrutura de resultados em três partes:** (1) validação da ferramenta, (2) aplicação às mutações de doença, (3) comparação com um controle neutro.
- **Considerar a especificidade de ligação:** eles também avaliaram mudanças em *qual* sequência de DNA a proteína reconhece, por meio de logos de ligação. Nós não fizemos isso, e isso entra como limitação e trabalho futuro.
- **Discutir as implicações terapêuticas:** mutações que só desestabilizam poderiam, em tese, ser compensadas por moléculas estabilizadoras. Mutações de contato, não.

---

## Referências citadas (conferir todas antes da submissão)

- Abramson J, et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*.
- Alibés A, Nadra AD, De Masi F, Bulyk ML, Serrano L, Stricher F. (2010). Using protein design algorithms to understand the molecular basis of disease caused by protein–DNA interactions: the Pax6 example. *Nucleic Acids Res* 38(21):7422–7431. doi:10.1093/nar/gkq683. (PDF em `dados_brutos/literatura/alibes2010_nar_pax6.pdf`, repositório DASH de Harvard, licença CC BY-NC.)
- Avsec Ž, et al. (2026). AlphaGenome. **[conferir título e revista]**
- Barrera LA, et al. (2016). Survey of variation in human transcription factors reveals prevalent DNA binding changes. *Science*.
- Cheng J, et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. **[conferir]**
- Dieckhaus H, et al. (2024). ThermoMPNN. *PNAS*. **[conferir]**
- Hovde S, et al. (2001). Estrutura do homeodomínio Msx-1 com DNA (PDB 1IG7). **[conferir]**
- Intarak N, Tongchairati K, Termteerapornpimol K, Chantarangsu S, Porntaveetus T. (2023). Tooth agenesis patterns and variants in PAX9: a systematic review. *Jpn Dent Sci Rev* 59:129–137. doi:10.1016/j.jdsr.2023.04.001.
- Jaganathan K, et al. (2019). Predicting splicing from primary sequence with deep learning (SpliceAI). *Cell*. **[conferir]**
- Liang J, Von den Hoff J, Lange J, Ren Y, Bian Z, Carels CEL. (2016). MSX1 mutations and associated disease phenotypes: genotype-phenotype relations. *Eur J Hum Genet* 24(12):1663–1670. doi:10.1038/ejhg.2016.78.
- Pak MA, et al. (2023). Using AlphaFold to predict the impact of single mutations on protein stability and function. *PLoS One*. **[conferir]**
- Pejaver V, et al. (2022). Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3 and BP4 criteria. *Am J Hum Genet*. **[conferir]**
- Peters H, et al. (1998). Pax9-deficient mice. *Genes Dev*. **[conferir]**
- Satokata I, Maas R. (1994). Msx1 deficient mice. *Nat Genet*. **[conferir]**
- Vastardis H, et al. (1996). A human MSX1 homeodomain missense mutation causes selective tooth agenesis. *Nat Genet*.
- Wang Y, et al. (2009). Pathogenic mechanisms of tooth agenesis linked to paired domain mutations in human PAX9. *Hum Mol Genet*. doi:10.1093/hmg/ddp221. PMC2706687.
- Xu HE, et al. (1999). Estrutura do domínio paired do PAX6 com DNA (PDB 6PAX). *Genes Dev*. **[conferir]**
- Revisão sobre a via Wnt e outras na agenesia não sindrômica (PubMed 29969831). **[conferir autores e ano]**
- Estudo de calibração PP3/BP4 com estrutura no BRCA1 (PMC12120176). **[conferir autores]**
- Estudo de variantes de homeodomínios (*Nature Communications*, 2024). **[conferir autores]**
