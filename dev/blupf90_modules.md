# BLUPf90 modules in `inst/software/` — what each does, when to use it, and how they connect

A reference for the executables bundled in `inst/software/`, grounded in the BLUPf90
manual (`inst/software/readme/blupf90_all8.pdf`, Misztal et al. 2022). Programs are
grouped by the stage of a genetic evaluation. "BIGf90" notes which existing package
function already drives that module (blank = no wrapper yet).

**One-line pipeline:**
raw data + pedigree + genotypes → *(QC)* → **renumf90** → *(preGSf90 for genomic)* →
**blupf90+** / **gibbsf90+ → postgibbsf90** → **postGSf90 / predf90 / predictf90 / validationf90**.

Two ground rules from the manual:

- The *application* programs (blupf90+, gibbsf90+, postgibbsf90, preGSf90, postGSf90,
  predictf90) are driven by a **parameter file** and need data with effects renumbered
  from 1. **renumf90** produces that renumbered set, so it runs first.
- A few tools (predf90, qcf90, seekparentf90) take **command-line arguments** instead
  of a parameter file.

---

## Stage 1 — Preparation & quality control (before renumbering)

| Module | What it does | Inputs | Outputs | When to use |
|---|---|---|---|---|
| **qcf90** | Stand-alone quality control on genotypes + pedigree (call rate, MAF, monomorphic, Mendelian, HWE). Command-line. | genotype file, pedigree, map | cleaned genotypes, QC reports | Clean SNP/pedigree data up front, independent of building G. |
| **seekparentf90** | Verifies paternity/maternity and discovers parents from SNP markers. Command-line. | SNP file, pedigree, map | parentage conflict + candidate-parent reports | Check/repair the pedigree with genotypes before evaluation. |
| **inbupgf90** | Inbreeding coefficients with incomplete pedigree; can assign unknown parent groups (UPG). | pedigree file | inbreeding (F) per animal, UPG assignments | Need F, or UPG for missing parents, before/around renumbering. |

*R side (BIGf90):* `write_par()` builds the raw `.par`; `write_geno()` builds the `.geno`.

---

## Stage 2 — Renumbering (always first for the application programs)

| Module | What it does | Inputs | Outputs | BIGf90 |
|---|---|---|---|---|
| **renumf90** | Recodes/renumbers effects, data and pedigree from 1; assigns UPG; basic pedigree/genotype checks. Writes the parameter + data files the application programs consume. | raw parameter file (your `.par`), data file, pedigree file, *(optional)* SNP file, *(optional)* map | `renf90.par` (renumbered par), `renf90.dat` (renumbered data), `renadd0X.ped` (renumbered pedigree + F for effect X), `renf90.tables` (recode keys), `<geno>_XrefID` (renumbered↔original genotype IDs), `renf90.inb` | `run_renum()` |

**When:** always, before any application program — everything downstream reads `renf90.par`
and the `renadd`/`_XrefID` files.

---

## Stage 3 — Genomic relationships (genomic models only)

| Module | What it does | Inputs | Outputs | When to use |
|---|---|---|---|---|
| **preGSf90** | Genomic preprocessor: QC of SNPs, builds and inverts **G** (genomic) and **A22** (pedigree for genotyped animals), and forms `G⁻¹ − A22⁻¹` so BLUP becomes single-step GBLUP (ssGBLUP). | `renf90.par` with `OPTION SNP_file`, renumbered genotypes (`_XrefID`), *(optional)* `map_file` | `GimA22i` (binary `G⁻¹−A22⁻¹`), `freqdata.count[.after.clean]` (allele freqs + exclusion codes), `Gen_call_rate`, `Gen_conflicts` (Mendelian conflicts), *(optional)* `G`/`A22`/inverses | Any GBLUP/ssGBLUP run. Usually invoked **automatically** by the solver when the par has `SNP_file`; run standalone to pre-build/save G or its QC. |

Key `OPTION`s: `whichG` (how G is formed), `whichfreq`/`FreqFile` (allele freqs),
`weightedG <file>` (weighted G = ZDZ′, weights from postGSf90), and QC knobs
`minfreq`, `callrate`, `callrateAnim`, `monomorphic`, `hwe`, `verify_parentage`.

---

## Stage 4 — Estimation: variance components & breeding values

| Module | What it does | Inputs | Outputs | When to use | BIGf90 |
|---|---|---|---|---|---|
| **blupf90+** | Combined blupf90 / remlf90 / airemlf90. Solves the mixed model for BLUP/GBLUP/ssGBLUP **solutions (EBVs)** and/or estimates **variance components** by EM-REML or AI-REML. | `renf90.par` *(+ `SNP_file` for genomic)* | `solutions` (trait effect level solution [SE]), variance-component estimates, *(optional)* `yhat_residual` | Get EBVs/GEBVs, and/or REML variances. AI-REML = fast + gives SE; EM-REML = robust/slow. `maxrounds 0` in AIREML → BLUP only. | `run_blup()`, `clean_ebvs()` |
| **gibbsf90+** | Combined Gibbs samplers (linear + threshold/categorical, homo/heterogeneous residuals). Bayesian variance components via MCMC. | `renf90.par` *(+ `SNP_file`)* | posterior sample files (`gibbs_samples`, `fort.99`…), `last_solutions` | Variance components by MCMC — complex, categorical, or over-parameterized models; cross-check on REML. | `run_gibbs()` |
| **postgibbsf90** | Post-processes the Gibbs chain: posterior means/SD and convergence diagnostics. | `gibbs_samples` + the parameter file | `postmean`, `postsd`, posterior summaries, convergence plots | Always after gibbsf90+, to summarize posteriors and check mixing/convergence. | `run_postgibbs()` |

---

## Stage 5 — Post-processing, prediction & validation

| Module | What it does | Inputs | Outputs | When to use | BIGf90 |
|---|---|---|---|---|---|
| **postGSf90** | Genomic postprocessor: back-solves **SNP effects** from GEBVs (ssGWAS), variance explained by SNP windows, SNP weights, and (optionally) p-values. | `renf90.par` + `solutions` (from a genomic blupf90+ run) + genotypes + `map_file` | `snp_sol` (SNP solution + weight + window variance), `chrsnp`/`chrsnpvar`, `windows_variance`, `snp_pred` (freqs + SNP effects), Manhattan plots, `snp_var` (with `snp_p_value`) | GWAS / SNP effects, variance by genomic window, or SNP weights to build a weighted G (iterate with preGSf90). | — |
| **predf90** | Predicts **direct genomic values (DGV = Zâ)** for young/new genotyped animals from existing SNP effects — no re-run. Command-line. | `snp_pred` (from postGSf90) + genotype file of new animals (`--snpfile`); `snp_var` for `--acc` | DGV predictions, *(optional)* reliabilities | Indirect prediction of animals genotyped after the main evaluation. | — |
| **predictf90** | Computes **y\*** (phenotype adjusted for the non-included effects), **ŷ** (included effects) and residuals, so `cor(ŷ, y*) ≈ accuracy`. | a blupf90 parameter file + `solutions` + data file + `OPTION include_effects` | `yhat_residual` (y\*, ŷ, residual), plus an EBV file if animal is in the model | Cross-validation / model accuracy — this is the engine BIGf90's `run_cva()` uses. | `run_predict()` |
| **validationf90** | Validation of predictions via the LR (Legarra–Reverter) method — bias, dispersion, and accuracy from comparing partial- vs whole-data solutions. | solutions from reduced vs full data (+ pedigree/genomic) | LR validation statistics (bias, slope/dispersion, accuracy ratio) | Forward/predictive validation across a data cut. *(Newer program; light coverage in the bundled manual.)* | — |
| **idsolf90** | Utility to attach **original IDs** to renumbered solutions (via `renadd0X.ped` / `_XrefID`). *(Not in the main manual's program list — described from its role.)* | `solutions` + `renadd0X.ped` / `_XrefID` | solutions keyed by original ID | Translate results back to real IDs. `clean_ebvs()` does the equivalent in R. | (`clean_ebvs()`) |

---

## Bundled helper scripts (not BLUPf90 executables)

| File | What it is |
|---|---|
| `blupf90_cv.par` | Example parameter file for the cross-validation workflow. |
| `cva_blupf90.sh` | Bash driver orchestrating k-fold CV (renumber → mask folds → blupf90+ → predictf90 → collect). `run_cva()` is the R reimplementation. |
| `acc_and_b1_linux.R` | Computes predictive accuracy and b1 (regression slope / bias) from the CV outputs. |
| `mcmc_coda_evaluation.R` | MCMC convergence diagnostics (via the `coda` package) for gibbsf90+ output. |

---

## How the modules connect

1. **Prepare inputs.** Build the model parameter file (`write_par()`) and, for genomic
   models, the genotype file (`write_geno()`). Optionally QC and fix the pedigree/genotypes
   first with **qcf90**, **seekparentf90**, and **inbupgf90**.
2. **Renumber.** **renumf90** turns your raw `.par` + data + pedigree (+ genotypes) into
   `renf90.par`, renumbered data/pedigree, and the `_XrefID` genotype key. Everything after
   this reads those files.
3. **Genomic setup (if applicable).** **preGSf90** QCs SNPs and forms `G⁻¹−A22⁻¹`
   (`GimA22i`) so the solver runs ssGBLUP. Triggered by `SNP_file` in the par.
4. **Estimate.** Run **blupf90+** for EBVs/GEBVs and AI/EM-REML variances, **or**
   **gibbsf90+** → **postgibbsf90** for Bayesian variances. Both produce a `solutions` file.
5. **Post-process / predict / validate.**
   - **postGSf90** back-solves SNP effects (ssGWAS) and SNP weights; those weights can feed
     **back** into preGSf90 as a weighted G (an iterative loop).
   - **predf90** turns SNP effects (`snp_pred`) into DGVs for newly genotyped animals.
   - **predictf90** produces y\* for cross-validation accuracy (used by `run_cva()`).
   - **validationf90** runs LR-method forward validation.
6. **Back to original IDs.** `clean_ebvs()` (or **idsolf90**) maps renumbered solutions to
   real animal IDs.

*Feedback loop:* SNP weights from **postGSf90** → weighted **G** in **preGSf90** → re-run
the solver (weighted ssGBLUP / iterative BLUP|SNP).
