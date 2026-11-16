# Agenesia dentária: levantamento de variantes PAX9 e MSX1

Atualizado em 15/11/2026.

## Arquivos

- `variantes_PAX9_MSX1.csv`: tabela consolidada, uma linha por troca de aminoácido, com ClinVar, gnomAD v4 e AlphaMissense.
- `dados_brutos/`: dados como vieram das fontes (ClinVar via NCBI E-utilities, gnomAD v4 via API, AlphaMissense via AlphaFold DB).
- `scripts/baixar_dados.py` e `scripts/consolidar.py`: refazem tudo do zero.

```
python -I scripts/baixar_dados.py dados_brutos
python -I scripts/consolidar.py dados_brutos variantes_PAX9_MSX1.csv
```

## Rótulos

- **P**: Pathogenic ou Likely pathogenic no ClinVar, sem conflito.
- **B**: Benign ou Likely benign no ClinVar, ou frequência no gnomAD de pelo menos 0,1%.
- **U**: o resto (incerto, conflitante, ou só raro no gnomAD).

Domínios de ligação ao DNA (UniProt): PAX9 paired 4–130; MSX1 homeobox 172–231.

## Primeiros números

| Gene | P | P no domínio | B | Incertas (U) | U no domínio |
|---|---|---|---|---|---|
| PAX9 | 19 | 16 | 4 | 531 | 147 |
| MSX1 | 5 | 2 | 7 | 775 | 72 |

AlphaMissense separa P de B com AUC 0,936 (24 P, 11 B).

## Ressalvas (importantes para o artigo)

1. **Confusão com o domínio.** No PAX9, todas as patogênicas com troca de aminoácido estão no domínio paired, e todas as benignas estão fora dele. Só "estar no domínio" já separa os grupos. O teste de verdade é dentro do domínio: as 147 + 72 incertas.
2. **Poucas benignas** (11), e algumas só por frequência no gnomAD.
3. **Circularidade.** Classificações recentes do ClinVar podem ter usado preditores computacionais (critério PP3), o que infla o AUC.
4. **Códon de início.** M1V, M1L e M1K do PAX9 são perda do códon de início, não trocas comuns: o AlphaMissense não foi feito para isso. Excluir da modelagem estrutural.
5. **Exceções no MSX1.** G122E e P247S são patogênicas no ClinVar, mas o AlphaMissense dá escore baixo. Vale olhar a literatura original.
6. **Revisão fraca.** Várias patogênicas estão com "no assertion criteria provided" ou têm um só submissor.

## Literatura (conferida em 15/11/2026)

Arquivo `literatura_missense.csv`, com 39 variantes:
- **PAX9:** 32 variantes, do Suplemento da revisão sistemática de Intarak et al. 2023 (Jpn Dent Sci Rev 59:129–137, doi:10.1016/j.jdsr.2023.04.001). Baixado pelo Europe PMC.
- **MSX1:** 7 variantes, da Tabela 2 de Liang et al. 2016 (Eur J Hum Genet 24:1663–1670, doi:10.1038/ejhg.2016.78).

Novo rótulo **P-lit**: variante descrita em família com agenesia, sem classificação patogênica no ClinVar, com frequência no gnomAD abaixo de 0,01% e sem conflito.

| Gene | P (ClinVar) | P-lit | P + P-lit no domínio |
|---|---|---|---|
| PAX9 | 19 | 19 | 32 |
| MSX1 | 5 | 6 | 7 |

Sem as perdas do códon de início, o AlphaMissense tem AUC 0,983 só com o ClinVar (21 P contra 11 B) e 0,977 com a literatura (44 contra 11). Valem as mesmas ressalvas acima.

Pontos a checar nos artigos originais:
- **PAX9 "Asn16Ile"** (Shahid 2016): na sequência de referência (P55771) a posição 16 é valina, não asparagina. É provável erro de numeração (talvez N17I) e ficou fora da tabela.
- **PAX9 Y143C**: conflitante no ClinVar, 0,035% no gnomAD, AlphaMissense 0,05. Provavelmente não é causal.
- **PAX9 A168G**: AlphaMissense 0,08. Os portadores tinham só 2 dentes ausentes. Evidência fraca.
- **MSX1 A225T**: recessiva (Chishti 2006).
- **MSX1 G122E e P247S**: patogênicas no ClinVar, mas com escore baixo. Ler os artigos.
- **LOVD:** o site bloqueou o IP em 15/11. Para usar, escrever para eles pedindo desbloqueio.

## AlphaFold Server (15/11/2026)

`jobs_afserver.csv` lista os 41 jobs (sequência, DNA e status). Gerado por `scripts/gerar_jobs.py`.

Construções:
- **PAX9:** resíduos 6–137 mais as duas fitas de DNA da estrutura 6PAX (PAX6).
- **MSX1:** resíduos 166–237 mais as duas fitas de DNA da estrutura 1IG7 (Msx-1 de camundongo).

Andamento:
- **Concluído:** PAX9_WT (no servidor, com o nome `PAX9_paired_domain_DNA_WT`), com ipTM 0,86 e pTM 0,85.
- **Enviados em 15/11 (29):** 21 mutantes do PAX9 (G6R até R77G), mais MSX1_WT e 7 mutantes do MSX1.
- **Pendentes para a cota de 16/11 (11):** PAX9 R77Q, T80A, I87F, K91E, A99P, A108P, V111M, C112W, V117G, S119F, S122Y.

Uma semente automática por job, com 5 amostras cada.

### Resultados dos 30 jobs (lote baixado em 15/11)

- O .zip em lote está em `dados_brutos/afserver/lote_parcial.zip` (57 MB).
- `scripts/extrair_afserver.py` lê o .zip sem extrair nada e escolhe, em cada job, a amostra de maior ranking_score.
- `scripts/distancia_dna.py` mede a distância de cada sítio mutado ao DNA na estrutura normal.
- Tudo vai para `resultados_afserver.csv`.

1. **A confiança global não discrimina.** Os ipTM vão de 0,84 a 0,87 e os pTM de 0,84 a 0,89, nos normais e nos mutantes, sem colisões entre átomos. O AlphaFold 3 prevê o mesmo dobramento com a mutação, como já se esperava. Não usar o ipTM como medida de efeito.
2. **O que informa é onde a mutação está:**
   - **PAX9**, 21 mutantes:
     - contato direto com o DNA (menos de 4 Å): 11, que são F15I, N17K, P20L, L21P, R26W, R44C, S49L, G51A/S/V e R77G;
     - próximos (4 a 8 Å): 6;
     - distantes (mais de 8 Å): 4, que são G6R, L27P, R28P e G64V. Provavelmente agem desestabilizando a proteína. Três deles põem prolina ou valina em regiões estruturadas.
   - **MSX1**, 7 mutantes:
     - contato: R202H e R202P;
     - próximos: 4;
     - distante: A227E.
   - Todos os sítios têm pLDDT entre 86 e 99, ou seja, estão em regiões bem ordenadas.
3. **Linha de base:** no domínio paired, 33% dos resíduos tocam o DNA. Entre os sítios patogênicos únicos do PAX9, 9 de 17 tocam (53%). Pelo teste binomial, p = 0,07, o que sugere enriquecimento mas ainda não é significativo. Faltam 11 mutantes do PAX9.

## AlphaGenome (15/11)

- **Instalação:** pacote `alphagenome` 2.0.0, num ambiente separado em `%USERPROFILE%\venvs\alphagenome`. A chave fica em `candidaturas-odontologia/.env`.
- **Variantes:** `scripts/clinvar_nao_codificantes.py` baixou do ClinVar 84 variantes fora da região codificante, com coordenadas no genoma de referência (GRCh38). O arquivo é `dados_brutos/clinvar_nao_codificantes.csv`.
- **Cálculo:** `scripts/alphagenome_score.py` usa os escores de splicing (sítios, uso de sítios e junções) e de expressão (RNA-seq), sempre do próprio gene. O AlphaGenome 2.0 exige escolher um tipo de célula, e não há tecido dentário. Usamos três: célula da crista neural (CL:0011012), célula-tronco mesenquimal (CL:0000134) e língua (UBERON:0001723). Fica o maior valor entre os três.
- **Piloto:** 7 patogênicas contra 7 benignas, em `alphagenome_piloto.csv`, com a saída bruta em `*_bruto_*.csv`. O escore de sítios de splicing separou os grupos por completo:
  - patogênicas: de 0,046 a 0,82, e seis delas acima de 0,6;
  - benignas: de 0,005 a 0,033.
  - A patogênica com escore baixo é a g.36662092G>A do PAX9, que provavelmente age por outro mecanismo. O RNA-seq dela (0,21) foi o mais alto do grupo.
  - Ressalvas: são poucas variantes, e as de splicing canônico (+1, +5) são casos fáceis.
- **Rodada completa:** saída em `alphagenome_todas.csv`. Das 81 variantes pontuadas, 78 receberam escore: 7 patogênicas, 37 benignas e 34 incertas. Três variantes intergênicas do MSX1 ficaram sem escore do próprio gene.
  - **Splicing:** AUC 0,988 entre patogênicas e benignas. A benigna mais alta chegou a 0,125, e a patogênica mais baixa ficou em 0,046.
  - **Expressão (RNA-seq):** AUC 0,985.
  - **Sete incertas com escore de splicing alto** (acima de 0,46, bem acima de todas as benignas) são candidatas a reclassificação:
    - MSX1 c.470-9G>A (0,93) e c.470-10C>A (0,92), duas trocas no íntron, perto do sítio aceptor;
    - PAX9 c.772-2A>G (0,76; sítio aceptor canônico) e c.772-3C>G (0,70);
    - PAX9 c.771+1G>A (0,64; sítio doador canônico), c.771+5G>C (0,61) e c.771+3G>T (0,46).
  - Ressalva: as variantes +1 e −2 são casos óbvios, que qualquer preditor de splicing pega. O artigo precisa comparar com um preditor padrão, como o SpliceAI. Os achados de verdade são os do MSX1 em −9 e −10.
  - Três variantes da região 3' não traduzida do PAX9 (c.*775, c.*231, c.*825) ficaram pouco acima da benigna mais alta. O sinal é fraco.

### AlphaGenome contra SpliceAI (15/11)

- **Método:** `scripts/spliceai_vep.py` busca os escores pré-calculados do SpliceAI pelo Ensembl VEP REST. O servidor do Broad não respondeu. Saída em `splicing_alphagenome_vs_spliceai.csv`.
- **Exclusões:** a PAX9 g.36662092G>A saiu da comparação. É perda do códon de início (M1I), não splicing. Duas deleções do MSX1 ficaram sem escore do SpliceAI.
- **Resultado com 5 patogênicas e 36 benignas:** os dois métodos têm AUC 1,000.
  - AlphaGenome: benigna mais alta 0,125; patogênica mais baixa 0,61.
  - SpliceAI: benigna mais alta 0,33; patogênica mais baixa 0,77.
- **Concordância nas incertas:** as 7 incertas que o AlphaGenome destacou também têm SpliceAI alto, entre 0,40 e 1,00. Os dois métodos dão 1,00 para MSX1 c.470-9G>A e c.470-10C>A, o que indica criação de um novo sítio aceptor. Abaixo disso, os dois concordam que o efeito é baixo, com uma exceção: na c.*775C>T o AlphaGenome dá 0,25 e o SpliceAI dá 0.
- **Conclusão honesta:**
  - Para splicing, o AlphaGenome confirma o SpliceAI, mas não acrescenta poder de discriminação neste conjunto.
  - Dois preditores independentes concordando reforçam o critério computacional (PP3) para essas 7 incertas, sobretudo as duas do MSX1.
  - O diferencial possível do AlphaGenome é o efeito na expressão de variantes reguladoras (região não traduzida e promotor). Isso o SpliceAI não cobre, mas ainda faltam patogênicas reguladoras confirmadas para validar.

### AlphaGenome 1 contra AlphaGenome 2 (15/11)

- **Método:** `scripts/alphagenome1_score.py` roda a API v1, que não pede tipo de célula e devolve escores normalizados (quantis). `scripts/comparar_versoes.py` faz a comparação, com saída em `comparacao_ag1_ag2_spliceai.csv`. Ficaram 6 patogênicas e 37 benignas, sem as perdas do códon de início.

| Método | AUC splicing | AUC expressão |
|---|---|---|
| AlphaGenome 2 (3 células) | 1,000 | 0,982 |
| AlphaGenome 1 (todas as células) | 1,000 | 0,959 |
| AlphaGenome 1 (as mesmas 3 células) | sem trilhas de splicing nessas células no v1 | 0,982 |
| SpliceAI (5 P, 36 B) | 1,000 | — |

- **Correlação de Spearman entre os métodos:** splicing v2 × v1, 0,70; v2 × SpliceAI, 0,50; v1 × SpliceAI, 0,58. Expressão v2 × v1, 0,82.
- **Ressalva de escala:** no v1, o escore em quantil satura perto de 1. As 7 incertas de splicing ficam todas em 1,00, e três da região 3' não traduzida chegam a 0,99. Essa escala é boa para ordenar as variantes, mas não para comparar valores absolutos com o v2 (escore bruto).
### Escolha da célula: mesenquimal como principal (15/11)

O tecido ideal seria o mesênquima do germe dentário, que não existe no AlphaGenome. Por isso a análise principal usa a **célula-tronco mesenquimal** (CL:0000134), que tem expressão e splicing medidos. A crista neural, a língua e o osteoblasto (UBERON:0002481_CL:0000062) entram como sensibilidade, cada uma relatada em separado. Scripts: `scripts/expressao_por_celula.py` e `scripts/analise_por_celula.py`. Saídas: `expressao_por_celula.csv` e `analise_por_celula.csv`.

**1. Os genes são expressos nessas células?** Usei o sinal de RNA previsto pelo AlphaGenome 2 no pico dentro do gene (percentil 99), dividido pelo sinal da albumina (ALB, que só aparece no fígado e serve de controle negativo).

| Gene | mesenquimal | crista neural | osteoblasto | língua |
|---|---|---|---|---|
| MSX1 | 35× | 19× | 55× | o controle ALB veio alto (1,95), então essa razão não serve |
| PAX9 | 11× | **0,7× (não expresso)** | 7× | sinal bruto 0,085 |

- Os dois genes têm expressão baixa em todas as células: o GAPDH tem sinal centenas de vezes maior. Era de esperar, porque nenhuma delas é germe dentário.
- **O PAX9 não aparece na crista neural.** Os escores do PAX9 nessa célula valem pouco.
- A **mesenquimal** é a única em que os dois genes ficam claramente acima do controle negativo. Isso confirma a escolha dela como principal.

**2. Desempenho por célula:** variantes patogênicas contra benignas, sem perda do códon de início. O limiar de cada célula é o maior escore de splicing entre as benignas.

| Célula | AUC splicing | AUC expressão | Limiar de splicing | Incertas acima do limiar |
|---|---|---|---|---|
| **mesenquimal (principal)** | 1,000 | 0,941 | 0,039 | **7** |
| osteoblasto | 1,000 | 0,928 | 0,054 | as mesmas 7 |
| língua | 1,000 | 0,986 | 0,125 | as 7, mais 3 da região 3' não traduzida (c.*825, c.*775, c.*231) |
| crista neural | 0,982 | 0,950 | 0,019 | 6, numa combinação diferente: perde as variantes do doador do PAX9 e acrescenta MSX1 c.470-5 e c.470-6 |

- **Conclusão:** as **7 incertas** (MSX1 c.470-9G>A e c.470-10C>A; PAX9 c.772-2A>G, c.772-3C>G, c.771+1G>A, c.771+5G>C e c.771+3G>T) se repetem na célula mesenquimal, no osteoblasto e na língua, e coincidem com o SpliceAI. São robustas.
- A crista neural diverge, o que bate com o PAX9 não ser expresso nela e com os dados de splicing dela serem estimados.
- As 3 variantes da região 3' não traduzida que só aparecem na língua têm evidência fraca.

- **Conclusão da comparação de versões:** as duas versões levam às mesmas 7 incertas de splicing. Em expressão, escolher células próximas ao dente melhora um pouco a discriminação (0,959 → 0,982), e o v2 faz isso de forma nativa. A troca de versão não muda nenhuma conclusão.

### Os 11 mutantes restantes do PAX9 (15/11, 21h30)

- **Envio:** R77Q, T80A, I87F, K91E, A99P, A108P, V111M, C112W, V117G, S119F e S122Y foram enviados depois que a cota renovou. Saiu um segundo lote em `dados_brutos/afserver/lote2.zip`, já incluído em `resultados_afserver.csv`. Os 41 jobs estão concluídos.
- **Confiança global:** ipTM de 0,84 a 0,87. A mesma ausência de discriminação de antes.
- **Posição em relação ao DNA:**
  - contato: A99P, R77Q, S119F, S122Y e V117G;
  - próximos: T80A e K91E;
  - distantes: I87F, A108P, V111M e C112W.
- **Total do PAX9, com 32 variantes em 27 sítios únicos:**
  - contato direto com o DNA: 13 sítios (48%), contra 33% no domínio todo. Teste binomial: p = 0,074;
  - até 8 Å do DNA: 19 sítios (70%), contra 66% no domínio. p = 0,40.
- **Conclusão honesta:** **o enriquecimento de contato com o DNA não é significativo** com 27 sítios. Quase metade das mutações patogênicas fica longe do DNA, no miolo da proteína, em geral trocas para prolina, metionina, triptofano ou fenilalanina em regiões estruturadas. Isso indica que o mecanismo predominante pode ser a **desestabilização da proteína**, não só a perda de contato com o DNA. A próxima análise deve estimar a estabilidade de cada mutante.

### Estabilidade (ΔΔG): tentativa com o DDMut (15/11, noite)

- `scripts/ddmut_estabilidade.py` extrai a cadeia proteica do modelo normal (`ddmut/PAX9_WT_cadeiaA.pdb`) e envia as mutações à API pública do DDMut.
- O envio em lista é recusado. O envio de uma mutação por vez é aceito, mas o resultado volta "Internal Server Error", **até no exemplo oficial da documentação (1BA3 H461D)**. O servidor está com falha, e não há ΔΔG calculado ainda.
- Opções: tentar o DDMut de novo mais tarde, ou rodar localmente uma ferramenta aberta, como o ThermoMPNN (Kuhlman Lab, licença MIT).

### Estabilidade (ΔΔG) com o ThermoMPNN (15/11, noite)

- **Instalação:** ThermoMPNN, clonado de github.com/Kuhlman-Lab/ThermoMPNN para `%USERPROFILE%\ferramentas\ThermoMPNN`, com os pesos `models/thermoMPNN_default.pt`. Roda em CPU, num ambiente separado (`%USERPROFILE%\venvs\thermompnn`, torch 2.5.1+cpu).
- **Entrada:** só a cadeia proteica do modelo normal do AlphaFold Server (`estruturas_wt/`). O programa calcula todas as trocas possíveis: 2.508 no PAX9 e 1.368 no MSX1. Saída em `estabilidade/`.
- **Análise:** `scripts/analise_estabilidade.py`. A convenção de sinal foi inferida pela mediana de todas as trocas (0,56 e 0,29, positiva), então ΔΔG maior que 0 significa desestabilizante.
- **Cuidado:** o ThermoMPNN numera as posições **a partir de 0**. Uma primeira rodada com a numeração errada foi descartada. O script agora confere o aminoácido de referência.

**Resultados do PAX9 (32 patogênicas):**
- A mediana do ΔΔG é 1,17, contra 0,56 entre todas as trocas possíveis (AUC 0,64).
- Por posição em relação ao DNA:

| Posição | n | Mediana do ΔΔG |
|---|---|---|
| contato com o DNA (menos de 4 Å) | 16 | **0,49**, igual ao fundo |
| próximas (4 a 8 Å) | 8 | **1,84** |
| distantes (mais de 8 Å) | 8 | **1,81** |

- As mutações **fora do contato desestabilizam mais** do que as de contato: AUC 0,88, teste de Mann-Whitney com p = **0,0003** (16 contra 16).
- **30 de 32** patogênicas se explicam por um dos dois mecanismos: contato com o DNA ou desestabilização (ΔΔG de pelo menos 1). Ficam sem explicação só G6R (0,44, distante) e R47W (0,62, a 4,5 Å).
- **Conclusão: são dois mecanismos complementares.** As mutações que tocam o DNA quase não desestabilizam, então devem atrapalhar a ligação. As que ficam longe do DNA desestabilizam o domínio. Isso fortalece a narrativa do artigo.
- **Ressalvas:**
  - o corte de 1 no ΔΔG foi escolhido depois de ver os dados;
  - o ThermoMPNN tem erro típico em torno de 1;
  - os valores vêm de uma única estrutura prevista, sem o DNA.

### Grupo controle do gnomAD (15/11, noite)

- **Método:** `scripts/grupo_controle.py`, com saída em `grupo_controle.csv`.
- **Controle:** variantes do gnomAD v4 dentro dos domínios (PAX9 6–130; MSX1 172–231) que trocam um aminoácido e **não** são patogênicas (nem P nem P-lit). Em duas versões: todas as observadas pelo menos uma vez, e só as recorrentes (vistas duas vezes ou mais).
- **Comparação:** distância ao DNA na estrutura normal, ΔΔG do ThermoMPNN e AlphaMissense. Contato medido pelo teste exato de Fisher; ΔΔG e AlphaMissense pelo teste de Mann-Whitney.

| | PAX9 patogênicas (32) | Controle, todas (124) | p | Controle, recorrentes (48) | p |
|---|---|---|---|---|---|
| Contato com o DNA | 50% | 26% | 0,011 | 17% | 0,0025 |
| ΔΔG mediano | 1,17 | 0,41 | 6×10⁻⁵ (AUC 0,73) | 0,34 | 4×10⁻⁵ (AUC 0,77) |
| Contato ou ΔΔG de pelo menos 1 | 94% | 40% | 1×10⁻⁸ | 33% | 3×10⁻⁸ |
| AlphaMissense (AUC) | — | 0,82 | 2×10⁻⁸ | 0,87 | 3×10⁻⁸ |

| | MSX1 patogênicas (7) | Controle, todas (59) | p | Controle, recorrentes (29) | p |
|---|---|---|---|---|---|
| Contato com o DNA | 29% | 24% | 1,0 | 21% | 0,64 |
| ΔΔG mediano | 1,52 | 0,16 | 0,004 (AUC 0,84) | 0,21 | 0,002 (AUC 0,88) |
| Contato ou ΔΔG de pelo menos 1 | 71% | 31% | 0,045 | 21% | 0,018 |

- **Conclusão:** frente às variantes da população, as patogênicas do PAX9 são enriquecidas **tanto em contato com o DNA quanto em desestabilização**. Os dois critérios juntos cobrem 94% delas contra 33–40% do controle. No MSX1, o sinal vem só da desestabilização: o contato não difere do controle, mas há só 7 patogênicas.
- **Ressalvas:**
  - o corte de 1 no ΔΔG foi definido depois de ver os dados. O ΔΔG contínuo, sozinho, já é significativo;
  - o controle pode conter variantes patogênicas raras ainda não descritas. Isso puxa o resultado para a hipótese nula, então o teste é conservador;
  - o AlphaMissense dá escore alto até para o controle (mediana de 0,94 a 0,99). O domínio inteiro é muito conservado, então o AlphaMissense discrimina menos aqui do que a combinação de estrutura e estabilidade.

### Melhoria 1: validação com experimentos publicados (15/11)

- **Fonte:** Wang et al. 2009, *Hum Mol Genet*, doi:10.1093/hmg/ddp221, PMC2706687. O estudo testou 8 mutantes do domínio paired: ligação ao DNA (EMSA), localização na célula, interação com o MSX1 e ativação dos promotores de Bmp4 e Msx1.
- **Arquivos:** `validacao_wang2009.csv` e `scripts/validacao_experimental.py`.
- ⚠️ **Os valores foram extraídos da página do PMC por leitura automática.** É preciso conferir na Tabela 1 e na Figura 3 do artigo original antes de publicar. O próprio artigo é inconsistente sobre a R26W: o texto diz que ela não liga ("−"), mas a legenda da Figura 3 diz "−/+".

**Regra de previsão:** a mutação perde a ligação ao DNA se fica em contato com o DNA ou tem ΔΔG de pelo menos 1.
- **Acertos: 6 de 8.** L21P, R26W, R28P, S43K e I87F perdem a ligação, como previsto. A G6R mantém a ligação, como previsto.
- **Os 2 erros são exceções conhecidas:**
  - a **G51S** liga normalmente ao DNA, mas perde a cooperação (sinergia) com o MSX1 e atrapalha a proteína normal (efeito dominante-negativo);
  - a **K91E** liga ao DNA, mas não ativa os promotores.

  As duas agem por mecanismos que não são ligação nem estabilidade, e nosso modelo não os cobre. Isso é uma limitação a declarar.
- **G6R:** foi a única sem nenhum dos dois mecanismos no nosso modelo. Ela também é a única que liga normalmente sem nenhum defeito descrito, e tem o quadro mais leve (3 dentes ausentes).
- **Correlações de Spearman (n = 8, valores apenas descritivos):**
  - ΔΔG do ThermoMPNN × estabilidade modelada por Wang: 0,55;
  - ligação ao DNA × dentes ausentes: −0,72 (menos ligação, mais dentes ausentes);
  - ΔΔG × dentes ausentes: 0,48.

### Melhoria 2: mecanismo × gravidade clínica (15/11)

- **Método:** `scripts/gravidade_clinica.py`. Usa o número de dentes ausentes, sem contar os terceiros molares, que vem da linha de cada variante no suplemento de Intarak 2023. Cruzei 25 mutações do PAX9 que trocam um aminoácido com o nosso mecanismo previsto. As 25 variantes nulas (frameshift e nonsense) entram como referência.

| Grupo | Dentes ausentes (mediana) | n |
|---|---|---|
| Contato com o DNA | 9,5 | 10 |
| Desestabiliza (ΔΔG de pelo menos 1, sem contato) | 10,0 | 13 |
| Nenhum dos dois | 9,5 | 2 |
| **Variantes nulas** | **15,0** | 25 |

- **O mecanismo não prevê a gravidade.** Desestabilizantes contra contato: AUC 0,57, p = 0,60. A correlação de Spearman entre ΔΔG e dentes ausentes é 0,23 (n = 25).
- **As mutações nulas são mais graves que as que trocam um aminoácido:** 15 contra cerca de 10 dentes, AUC 0,77, **p = 0,001**. Isso sugere que as trocas de aminoácido, pelos dois mecanismos, deixam **alguma função residual** (hipomórficas). As nulas eliminam a função de uma cópia inteira do gene.
- **Ressalvas:**
  - é um único valor por variante, de famílias com critérios de diagnóstico diferentes;
  - a mesma variante pode variar entre portadores;
  - a G64V, com 22 dentes ausentes, e a V117G, com 2, são extremos que pesam muito numa amostra pequena.

### Melhoria 3: interface com o DNA nos modelos mutantes (15/11)

- **Método:** `scripts/interface_dna.py`. Conta os contatos entre proteína e DNA (pares de átomos a menos de 4 Å) na melhor amostra de cada job e compara com a proteína normal. Também usa a diferença de ipTM entre proteína e DNA. Os resultados estão em `resultados_afserver.csv`, colunas `contatos_dna`, `delta_contatos_pct` e `delta_iptm_dna`.
- **Valores de referência da proteína normal:** PAX9 com 361 pares de contato e 44 resíduos na interface; MSX1 com 179 pares e 22 resíduos.
- **Resultado negativo:**
  - os contatos mudam entre −6% e +9% no PAX9, quase sempre para mais;
  - o ipTM entre proteína e DNA varia no máximo 0,035.
  - Mutantes que **perdem a ligação ao DNA no experimento** de Wang 2009 aparecem iguais ou com mais contatos: L21P +9,1%, R26W +0,3%, R28P +0,6%, I87F +5,8%.
- **Conclusão:** o AlphaFold 3 **não reproduz a perda de ligação causada por uma única mutação**. Ele "força" o complexo a se formar, porque aprendeu a partir de complexos parecidos. Isso confirma a limitação descrita na literatura e é um achado metodológico útil para o artigo. O efeito na ligação precisa vir da posição do resíduo (contato) e do ΔΔG, não da comparação entre modelos.

### Melhoria 4: classificação das variantes incertas dos domínios (15/11)

- **Método:** `scripts/classificar_incertas.py`, com saída em `incertas_classificadas.csv`.
- **Regra:** contato com o DNA ou ΔΔG de pelo menos 1. Calibrada com as patogênicas contra o controle do gnomAD (variantes recorrentes, vistas duas ou mais vezes). Para cada resultado da regra, calculei a razão de verossimilhança (LR): quanto o resultado aumenta (LR+) ou diminui (LR−) a chance de a variante ser patogênica.

| Gene | Sensibilidade | Falso-positivo | LR+ | LR− |
|---|---|---|---|---|
| PAX9 | 0,94 | 0,33 | 2,8 | **0,09** |
| MSX1 | 0,71 | 0,21 | 3,5 | 0,36 |

- **Leitura:** no PAX9, a regra é mais útil quando dá **negativa**. Com LR− de 0,09, ela é uma evidência forte de que a variante tende a ser benigna. Quando dá positiva, com LR+ de 2,8, a evidência é só de apoio, como o critério PP3 do ACMG.
- **Incertas nos domínios:**
  - PAX9: 138, sendo 59 com regra positiva e 79 com regra negativa;
  - MSX1: 72, sendo 24 positivas e 48 negativas.
- **Destaques entre as incertas registradas no ClinVar:**
  - **PAX9 positivas:** trocas para prolina ou arginina em regiões estruturadas, como L32P (ΔΔG 3,08), L106R (3,19), V83G (3,04), R105P e R88P. Também a R26Q, na mesma posição da patogênica R26W.
  - **MSX1 positivas:** W219R (conflitante no ClinVar; ΔΔG 2,61, encostada no DNA), R224C, A227P e K228E.
  - **PAX9 negativas:** D94E, A108V e Y115H, por exemplo. O AlphaMissense discorda em várias delas, como a C112F e a E62K, com escore de 1,00.
- **Ressalvas:**
  - a calibração foi feita na mesma amostra em que a regra é aplicada, então é otimista;
  - os valores de LR são aproximados;
  - não é diagnóstico clínico. Os termos do AlphaGenome e o bom senso pedem uso só em pesquisa.

### Melhoria 5: robustez da regra (15/11)

- **Método:** `scripts/robustez.py`. Só o PAX9: 32 patogênicas contra 48 variantes recorrentes do gnomAD.
- **(a) Validação cruzada deixando uma de fora (leave-one-out) para o corte de ΔΔG:** em todas as rodadas o corte escolhido foi **1,5**. Fora da amostra:
  - sensibilidade **0,84** e especificidade **0,79**;
  - **LR+ 4,0** e **LR− 0,20**.

  É a estimativa honesta para o artigo, porque é menos otimista que a calibração na própria amostra (0,94 e 0,67).
- **(b) Sensibilidade aos cortes**, sensibilidade / especificidade:

| Contato | ΔΔG ≥ 0,5 | ΔΔG ≥ 1,0 | ΔΔG ≥ 1,5 |
|---|---|---|---|
| < 3,5 Å | 0,88 / 0,62 | 0,81 / 0,71 | 0,72 / 0,83 |
| < 4,0 Å | 0,97 / 0,58 | 0,94 / 0,67 | 0,84 / 0,79 |
| < 5,0 Å | 0,97 / 0,54 | 0,97 / 0,62 | 0,88 / 0,75 |

  O resultado é estável: nenhuma combinação de cortes muda a conclusão.
- **(c) As 5 amostras do modelo normal:** 12 dos 13 sítios patogênicos em contato com o DNA aparecem em contato nas 5 amostras. Só a posição 117 (V117G) varia entre amostras. O contato é robusto.
- **Recomendação para o artigo:** usar o corte de **ΔΔG ≥ 1,5**, escolhido pela validação cruzada, com contato a menos de 4 Å. Reportar os valores fora da amostra (LR+ 4,0, LR− 0,20), e não os da calibração na própria amostra.

**Resultados do MSX1 (7 patogênicas):** a mediana do ΔΔG é 1,52, contra 0,29 entre todas as trocas (AUC 0,76). L211R, T180I, L230P e R202P estão entre as 16% trocas mais desestabilizantes do domínio. A R202H está em contato com o DNA e desestabiliza pouco (0,80). Com só 7 variantes, não dá para testar estatisticamente.

Próximo passo analítico: comparar com um grupo controle de variantes do gnomAD dentro do domínio (frequentes ou benignas) e estimar a estabilidade de cada mutante.

## Próximos passos

- AlphaFold Server: as 39 patogênicas (P + P-lit) dentro dos domínios (32 PAX9, 7 MSX1), mais a proteína normal do MSX1 ligada ao DNA. São uns 40 jobs, dois dias de cota. Depois, um grupo das incertas no domínio.

## 16/11/2026 — VUS reclassificadas com o limiar validado (ΔΔG ≥ 1,5)

`scripts/classificar_incertas.py` agora aceita o limiar como 6º argumento (padrão 1,5). Regra: contato < 4 Å ou ΔΔG ≥ 1,5.
- PAX9: 138 incertas no domínio, 45 positivas (LR+ 4,0) e 93 negativas (LR− 0,20). Antes, com ΔΔG ≥ 1: 59 positivas.
- MSX1: 72 incertas, 23 positivas (LR+ 3,5) e 49 negativas (LR− 0,36). Antes: 24 positivas.
- Exemplos do texto continuam positivos: L32P, L106R, V83G, R26Q, W219R, R224C, A227P.
- Versão antiga salva em `incertas_classificadas_ddg1.csv`.

## 16/11/2026 — auditoria para a passagem ao Codex

Ver `PASSAGEM_METODOLOGICA.md`, que **prevalece** sobre as seções antigas deste diário.
- Wang 2009 com ΔΔG ≥ 1,5: **5/8**. A S43K, com ΔΔG 1,46, passa a errar. O 6/8 acima vale só para ΔΔG ≥ 1.
- Conjunto patogênico do PAX9: 16 P (ClinVar) + 16 P-lit.
- Controles do gnomAD estão todos dentro das VUS reclassificadas (PAX9 124/124; recorrentes 48/48). Por isso a calibração e a aplicação não são independentes.
- IC de Wilson: sensibilidade 0,68–0,93; especificidade 0,66–0,88.
