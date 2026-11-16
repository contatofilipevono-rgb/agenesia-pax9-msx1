# MANUSCRIPT SKELETON — target: European Journal of Oral Sciences (plan B: Archives of Oral Biology)

*Draft 15 Nov 2026; corrected 16 Nov. **`PASSAGEM_METODOLOGICA.md` takes precedence over this skeleton wherever they differ** (validation, independence of sets, limitations). Numbers come from `LEIAME.md` and the CSV files in `agenesia/`. Items marked **[TODO]** or **[check]** must be resolved before submission. Check the journal's author guidelines for word limits, abstract format and reference style.*

---

## Title (options)

1. **Two complementary structural mechanisms explain most PAX9 missense variants in non-syndromic tooth agenesis: an AlphaFold 3, stability and AlphaGenome analysis**
2. Loss of DNA contact or domain destabilization: a structure-based framework to interpret PAX9 and MSX1 variants in tooth agenesis
3. Structural and regulatory interpretation of PAX9 and MSX1 variants in tooth agenesis using AlphaFold 3 and AlphaGenome

**Authors:** Filipe Martins Vono (Centro Universitário de Belo Horizonte – UniBH, School of Dentistry) **[TODO: co-authors — genetics/orthodontics; ORCID; corresponding author]**

---

## Abstract (structured, ~250 words — [check limit])

**Background.** Missense variants in the DNA-binding domains of PAX9 and MSX1 cause non-syndromic tooth agenesis, but how they disrupt function, and how to interpret new variants of uncertain significance (VUS), remains unclear.

**Methods.** We compiled 39 pathogenic missense variants in the PAX9 paired domain (n = 32) and the MSX1 homeodomain (n = 7) from ClinVar and the literature. We modelled wild-type and mutant protein–DNA complexes with AlphaFold 3 (41 predictions), estimated folding stability changes (ΔΔG) with ThermoMPNN, and compared pathogenic variants with population variants from gnomAD v4 in the same domains. Non-coding variants (n = 78) were scored with AlphaGenome 1 and 2 and compared with SpliceAI.

**Results.** AlphaFold 3 confidence did not change with mutation (ipTM 0.84–0.87). Pathogenic PAX9 variants were enriched both at DNA-contacting residues (50% vs 17%; p = 0.0025) and among destabilizing changes (median ΔΔG 1.17 vs 0.34 kcal/mol; p = 4×10⁻⁵). DNA-contact and destabilizing variants were mutually exclusive in mechanism (ΔΔG of non-contact vs contact variants, AUC 0.88; p = 0.0003). A combined rule (DNA contact or ΔΔG ≥ 1.5) achieved out-of-sample sensitivity of 0.84 and specificity of 0.79, and agreed with published DNA-binding assays for 5 of 8 mutants. The mechanism did not predict the number of missing teeth, but missense variants were milder than null variants (≈10 vs 15 teeth; p = 0.001). AlphaGenome and SpliceAI separated pathogenic from benign non-coding variants equally well (AUC 1.0) and jointly flagged seven splice-altering VUS, including two intronic MSX1 variants.

**Conclusions.** Most PAX9 agenesis variants act either by abolishing DNA contacts or by destabilizing the paired domain. A simple structural rule provides calibrated, research-level evidence for VUS interpretation.

**Keywords:** hypodontia; oligodontia; PAX9; MSX1; AlphaFold; protein stability; variant interpretation

---

## 1. Introduction (~600 words) — base: `revisao_literatura.md`

**¶1. Clinical problem.** Tooth agenesis: the most common dental developmental anomaly, with a strong genetic basis. ~200 variants in 15 genes, mostly AXIN2, EDA, LRP6, MSX1, PAX9, WNT10A and WNT10B [ref Wnt review **check**]. PAX9 is the most specific gene for non-syndromic forms [ref].

**¶2. Biology.** PAX9 and MSX1 act in dental mesenchyme at bud stage, sustain *Bmp4*, and arrest at bud stage when knocked out [Satokata & Maas 1994; Peters 1998 **check**]. They synergize on the *Bmp4* and *Msx1* promoters [Wang 2009].

**¶3. Genotype–phenotype and the gap.** Intarak et al. 2023 (61 PAX9 variants), Liang et al. 2016 (MSX1). Functional data exist for few variants (Wang 2009: 8 PAX9 mutants). Most missense variants in ClinVar are VUS (138 in the PAX9 paired domain, 72 in the MSX1 homeodomain in our data).

**¶4. Structural precedent.** In PAX6, Alibés et al. 2010 explained 88% of paired-domain missense variants by stability or DNA-binding changes (FoldX). No equivalent systematic analysis exists for tooth agenesis genes, and none uses current tools (AlphaFold 3, ThermoMPNN, AlphaGenome) or population controls.

**¶5. Aims.** (i) Classify PAX9/MSX1 pathogenic missense variants by mechanism (DNA contact vs destabilization) and test against population variants. (ii) Validate against published experiments and clinical severity. (iii) Assess AlphaGenome for non-coding variants. (iv) Provide calibrated evidence for VUS.

---

## 2. Material and methods (~1200 words)

### 2.1 Variant curation
- ClinVar (NCBI E-utilities, accessed 15 Nov 2026), PAX9 and MSX1, single-gene records. Missense matched to the canonical isoforms (UniProt P55771, P28360).
- Literature: Intarak et al. 2023 supplementary table (PAX9, 32 missense) and Liang et al. 2016 Table 2 (MSX1, 7). "P-lit" label: reported in an agenesis family, gnomAD AF < 10⁻⁴, not conflicting in ClinVar.
- gnomAD v4 (API, exomes and genomes); AlphaMissense (AlphaFold DB).
- Domains: PAX9 paired 4–130 (analysed 6–130); MSX1 homeobox 172–231.
- Exclusions: start-loss variants (M1); PAX9 "Asn16Ile" (inconsistent with the reference sequence).
- Final sets: 39 pathogenic (32 PAX9, 7 MSX1); controls = gnomAD missense in the domains, not pathogenic: all (PAX9 124, MSX1 59) and recurrent, AC ≥ 2 (PAX9 48, MSX1 29).

### 2.2 Structural modelling (AlphaFold 3)
- AlphaFold Server (15 Nov 2026), default settings, automatic seed, 5 samples per job; top-ranked sample used.
- PAX9 residues 6–137 with the double-stranded DNA from PDB 6PAX (PAX6 paired domain); MSX1 residues 166–237 with DNA from PDB 1IG7.
- 41 jobs (2 wild type, 39 mutants). Recorded ipTM, pTM, protein–DNA chain-pair ipTM, minimum PAE and clashes.
- Distance to DNA: minimum heavy-atom distance of each residue to DNA in the wild-type model; contact < 4 Å, near 4–8 Å, distant > 8 Å. Contact was stable across the 5 samples for 12/13 pathogenic contact sites.
- Interface in mutant models: protein–DNA atom pairs < 4 Å vs wild type.

### 2.3 Stability (ΔΔG)
- ThermoMPNN (Kuhlman Lab, default weights, CPU) run on the wild-type protein chain; site-saturation mutagenesis (2,508 PAX9 and 1,368 MSX1 substitutions). Sign convention checked (median of all substitutions > 0 = destabilizing). **[check: sign convention in the ThermoMPNN paper]**
- DDMut attempted but unavailable (server error, including the documentation example).

### 2.4 Statistics
- Fisher's exact test (contact proportion); Mann–Whitney U with AUC (ΔΔG, AlphaMissense); binomial test against the domain-wide contact fraction; Spearman correlation.
- Combined rule: DNA contact (< 4 Å) or ΔΔG ≥ threshold. Threshold chosen by leave-one-out cross-validation maximizing Youden's index (grid 0–3, step 0.25), then sensitivity, specificity and likelihood ratios out of sample. Sensitivity analysis: 3.5 / 4 / 5 Å and ΔΔG 0.5 / 1.0 / 1.5.
- Code in Python 3.12 (scripts in the repository **[TODO: GitHub/Zenodo DOI]**).

### 2.5 External validation
- Experimental: Wang et al. 2009 (EMSA, transactivation) for 8 PAX9 mutants. **[TODO: verify values against Table 1 and Fig. 3 of the original paper]**
- Clinical severity: number of missing teeth (excluding third molars) per variant from Intarak et al. 2023; null variants (frameshift/nonsense) as reference.

### 2.6 Non-coding variants (AlphaGenome)
- 84 non-coding ClinVar variants with GRCh38 coordinates (SPDI); 78 scorable after excluding protein-altering indels and empty alleles.
- AlphaGenome 2 (package `alphagenome` 2.0.0, API v2, default model, all folds), 1 Mb window, recommended scorers for splice sites, splice-site usage, splice junctions and RNA-seq, restricted to the variant's gene.
- Cell conditions (no dental tissue available): mesenchymal stem cell CL:0000134 (primary), neural crest cell CL:0011012, tongue UBERON:0001723 and osteoblast UBERON:0002481_CL:0000062 (sensitivity). Expression of PAX9/MSX1 in each condition checked with predicted RNA-seq against GAPDH (positive) and ALB (negative).
- AlphaGenome 1 (API v1, all tracks) for version comparison; SpliceAI precomputed scores via Ensembl VEP REST.

### 2.7 Ethics
- Public, de-identified data only; no human participants. **[check: whether the institution requires an ethics exemption statement]**

---

## 3. Results (~1200 words)

### 3.1 Dataset (Table 1)
- 39 pathogenic missense variants: 32 PAX9 (16 ClinVar P/LP + 16 literature-only, all in the paired domain) and 7 MSX1 (2 ClinVar + 5 literature-only). *(Corrected 16 Nov: the earlier "19 + 13" was wrong.)*

### 3.2 AlphaFold 3 confidence does not discriminate mutants (Fig. 1A, Suppl. Table)
- Wild type: PAX9 ipTM 0.86 / pTM 0.85; MSX1 0.87 / 0.88. All mutants: ipTM 0.84–0.87, no clashes.
- Mutant interfaces were unchanged (PAX9 contacts −6% to +9%). Experimentally DNA-binding-dead mutants (L21P, R26W, R28P, I87F) showed no interface loss (+0.3% to +9.1%).

### 3.3 Two complementary mechanisms (Fig. 2)
- **PAX9** (27 unique sites): 13 contact DNA (48% vs 33% domain-wide; binomial p = 0.074). By variant: 16 contact, 8 near, 8 distant.
- ΔΔG: contact median 0.49 (similar to background 0.56) vs near 1.84 and distant 1.81. Non-contact > contact: AUC 0.88, p = 0.0003.
- Examples. Destabilizing, non-contact: R28P 3.26, A108P 2.50, L27P 2.45, I56N 2.38. Contact, non-destabilizing: R26W, S49L, G51A/S/V, R77Q/G.
- **MSX1**: pathogenic median ΔΔG 1.52 vs 0.29 for all substitutions (AUC 0.76); L211R, T180I, L230P and R202P are in the top 16%.

### 3.4 Comparison with population variants (Table 2, Fig. 3)
- PAX9 vs recurrent gnomAD (AC ≥ 2, n = 48):
  - contact 50% vs 17% (p = 0.0025);
  - ΔΔG 1.17 vs 0.34 (AUC 0.77, p = 4×10⁻⁵);
  - contact or ΔΔG ≥ 1.5 (article rule): 84% vs 21% (Fig. 3B); with ΔΔG ≥ 1: 94% vs 33% (p = 3×10⁻⁸);
  - AlphaMissense AUC 0.87, with high scores even in controls (median 0.98).
- Same direction with all gnomAD variants (n = 124).
- MSX1: ΔΔG significant (p = 0.002); contact not different (2/7).

### 3.5 Cross-validated rule (Fig. 3C)
- Leave-one-out selected ΔΔG ≥ 1.5 in every fold. Out of sample: sensitivity 0.84, specificity 0.79, LR+ 4.0, LR− 0.20.
- Robust to distance and ΔΔG thresholds (Suppl. Table).

### 3.6 Experimental validation (Table 3)
- Wang et al. 2009, with the article rule (ΔΔG ≥ 1.5): 5/8 concordant. Correct: L21P, R26W, R28P, I87F (loss) and G6R (retained, mildest phenotype, 3 teeth).
- Discordant:
  - S43K: 4–8 Å from DNA, ΔΔG 1.46, just below the cut; reduced binding in EMSA.
  - G51S: binding retained, no MSX1 synergy, dominant-negative.
  - K91E: binding retained, no transactivation.
  - G51S and K91E act through mechanisms outside the model.
- With ΔΔG ≥ 1.0 the result is 6/8 (S43K concordant); report this as a sensitivity analysis.
- All 8 variants belong to the calibration set, so this validates the mechanism, not the classifier (see `PASSAGEM_METODOLOGICA.md` §4).

### 3.7 Clinical severity (Fig. 4)
- Missing teeth, ΔΔG ≥ 1.5 rule (Fig. 4): contact 9.5 (n = 10), destabilizing 11.0 (n = 11), neither 7.0 (n = 4); destabilizing vs contact p = 0.40. With ΔΔG ≥ 1: 9.5 vs 10.0, p = 0.60. Spearman ΔΔG vs teeth 0.23.
- Null variants 15.0 (n = 25) vs missense: AUC 0.77, p = 0.001.

### 3.8 Non-coding variants (Fig. 5)
- AlphaGenome 2 splice-site score: AUC 1.0 (6 pathogenic vs 37 benign); RNA-seq AUC 0.98. SpliceAI AUC 1.0. AlphaGenome 1 led to the same conclusions.
- Seven VUS flagged by all methods: MSX1 c.470-9G>A and c.470-10C>A (SpliceAI 1.00; likely new acceptor site); PAX9 c.772-2A>G, c.772-3C>G, c.771+1G>A, c.771+5G>C and c.771+3G>T.
- Cell choice: PAX9 not expressed in the neural crest condition (0.7× ALB); mesenchymal cells express both genes (PAX9 11×, MSX1 35× ALB); the 7 VUS replicate in mesenchymal, osteoblast and tongue.

### 3.9 VUS in the DNA-binding domains (Suppl. Table)
- Rule = DNA contact < 4 Å or ΔΔG ≥ 1.5 (cross-validated threshold; recomputed 16 Nov 2026, `incertas_classificadas.csv`).
- PAX9: 138 VUS, 45 rule-positive (LR+ 4.0) and 93 rule-negative (LR− 0.20).
- MSX1: 72 VUS, 23 rule-positive (LR+ 3.5) and 49 rule-negative (LR− 0.36; MSX1 calibration from only 7 pathogenic variants).
- Previous version with ΔΔG ≥ 1 (59 and 24 positive) kept in `incertas_classificadas_ddg1.csv`.
- Examples: PAX9 L32P, L106R, V83G, R26Q (same site as R26W); MSX1 W219R, R224C, A227P.

---

## 4. Discussion (~1200 words)

- **¶1 Main finding.** Two mechanisms explain most PAX9 variants, mirroring PAX6 (Alibés 2010: 88%; stability 77%, binding 30%). Our cross-validated threshold (1.5 kcal/mol) converges with their 1.6 kcal/mol (twice the FoldX error), and our population control is stronger than an ortholog set.
- **¶2 Mechanistic interpretation.** Contact variants cluster in the N-terminal subdomain β-hairpin and linker and spare stability. Distant variants introduce prolines or bulky residues into the helical cores (R28P, L27P, A108P). This is consistent with Wang 2009 (loss of binding without loss of nuclear localization) and with haploinsufficiency.
- **¶3 Exceptions show the model's limits.** G51S (dominant-negative, loss of MSX1 synergy) and K91E (transactivation). Interface specificity changes (Alibés' DNA logos) were not modelled. Homeodomain specificity variants are poorly predicted by current tools (Nat Commun 2024 **check**).
- **¶4 AlphaFold 3 is insensitive to point mutations.** Point mutants keep native-like interfaces (Pak 2023 **check**), so mutant models should not be used to infer loss of binding. Position and ΔΔG should be used instead.
- **¶5 Severity.** No mechanism–severity correlation (as in PAX6), but missense variants are milder than null variants, suggesting hypomorphic effects. Ascertainment and variability between carriers limit this analysis.
- **¶6 Non-coding variants and AlphaGenome.** AlphaGenome matches SpliceAI for splicing; its added value lies in regulatory effects and cell-type conditioning. There is no dental tissue among its conditions; future work could use dental pulp stem cell or tooth germ single-cell data with AlphaGenome 2 cell conditioning.
- **¶7 Clinical translation.** Likelihood ratios (LR+ 4.0, LR− 0.20) correspond to supporting or moderate evidence under the ClinGen PP3/BP4 framework (Pejaver 2022 **check**). They are for research use only, not diagnosis. Rule-negative results are the more informative (benign-supporting).
- **¶8 Limitations.**
  - Single predicted structure, without DNA in the ΔΔG calculation.
  - ThermoMPNN error ≈ 1 kcal/mol.
  - Small MSX1 sample.
  - Controls may hide rare pathogenic variants (conservative bias).
  - Experimental values taken from a single study.
  - Literature-derived severity data.
  - No wet-lab validation; splice VUS need minigene assays.
- **¶9 Conclusion.**

---

## 5. Figures and tables (plan)

| # | Content | Source |
|---|---|---|
| **Fig. 1** | (A) AlphaFold 3 models of PAX9–DNA and MSX1–DNA coloured by pLDDT; (B) pathogenic sites on the structure coloured by mechanism (contact = blue, destabilizing = red, neither = grey) | `dados_brutos/afserver/*.zip`, `resultados_afserver.csv` |
| **Fig. 2** | Scatter of distance to DNA (x) vs ΔΔG (y) for pathogenic PAX9 variants, with the two-mechanism quadrants | `resultados_afserver.csv` |
| **Fig. 3** | (A) ΔΔG distributions, pathogenic vs gnomAD; (B) contact proportions; (C) cross-validated rule performance | `grupo_controle.csv`, `scripts/robustez.py` |
| **Fig. 4** | Missing teeth by mechanism and null variants | `scripts/gravidade_clinica.py` |
| **Fig. 5** | AlphaGenome vs SpliceAI for non-coding variants; the 7 flagged VUS | `comparacao_ag1_ag2_spliceai.csv`, `analise_por_celula.csv` |
| **Table 1** | Pathogenic variants: source, distance, mechanism, ΔΔG, AlphaMissense | `resultados_afserver.csv` |
| **Table 2** | Pathogenic vs gnomAD statistics | `grupo_controle.csv` |
| **Table 3** | Concordance with Wang 2009 | `validacao_wang2009.csv` |
| **Suppl.** | Full VUS classification; AlphaGenome by cell type; sensitivity analysis; AlphaFold Server job list | CSV files in `agenesia/` |

**[TODO: generate figures (Python/matplotlib) and render the structures (PyMOL/ChimeraX, or Mol* screenshots)]**

---

## 6. Declarations

- **Data and code availability:** scripts and processed tables on GitHub/Zenodo **[TODO]**. AlphaFold Server outputs are shared under the AlphaFold Server Output Terms (non-commercial use, notice required). The AlphaGenome API key is not shared.
- **Tool terms:** AlphaFold Server and AlphaGenome are for non-commercial use only. AlphaGenome outputs must not be used for clinical purposes; state this in the Methods.
- **Funding:** none (or Google Cloud Research Credits, if awarded). **[TODO]**
- **Conflicts of interest:** none. **[check]**
- **Author contributions (CRediT):** **[TODO]**
- **Preprint:** bioRxiv or medRxiv before submission (recommended).

---

## 7. Pending before submission (checklist)

- [ ] Check the **[check]** references in `revisao_literatura.md` and here.
- [ ] Confirm the Wang 2009 values in Table 1 and Fig. 3 of the original paper.
- [x] Recompute the VUS classification with the cross-validated ΔΔG ≥ 1.5 threshold (16 Nov: PAX9 45/138, MSX1 23/72).
- [ ] Confirm the ThermoMPNN sign convention in its paper.
- [x] Generate figures 1–5 (`figuras/`).
- [ ] Define co-authors (genetics/orthodontics) and affiliations.
- [ ] Deposit code and data (GitHub + Zenodo DOI).
- [ ] Read the EJOS author guidelines (limits, abstract format, references).
- [ ] Optional: retry DDMut as a second ΔΔG predictor (robustness).
