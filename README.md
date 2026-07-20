# LN-PriGI-PRS

LN-PriGI-PRS is a standalone command-line tool for calculating a lupus-nephritis polygenic score (PGS) and reporting each sample's percentile relative to a specified reference population.

The percentile can support research on lupus-nephritis susceptibility. This software does **not** provide clinical risk interpretation, diagnosis, prognosis, or treatment advice.

## Requirements

- R 4.3 or later
- R packages: `data.table` and `digest`
- PLINK 1.9 installed locally and available on `PATH`, or supplied with `--plink`
- Genotypes aligned to the hg19/GRCh37 assembly

PLINK 2, PGEN/PVAR/PSAM, and VCF inputs are not supported. No liftover is performed.

## Input

Provide a file prefix for exactly one complete PLINK 1 file set:

- Binary: `PREFIX.bed`, `PREFIX.bim`, and `PREFIX.fam`
- Text: `PREFIX.ped` and `PREFIX.map`

The program accepts only biallelic variants. A/T and C/G palindromic SNPs are excluded. Variants are aligned by chromosome, hg19 position, and alleles; input variant IDs are not treated as cross-dataset identifiers.

Prepare and impute genotypes before running this tool. The program does not perform genotype imputation. A missing target genotype contributes zero to the raw PGS, and coverage is calculated from the genotype data before this scoring representation is created.

## Basic usage

```bash
Rscript R/main.R \
  --input /path/to/cohort \
  --genome-build hg19 \
  --output-dir /path/to/results \
  --output-prefix ln_prs \
  --plink /path/to/plink
```

The score is the uncentered, unstandardized sum produced by PLINK 1.9 `--score ... header sum`. `--maf` and `--geno` are optional conveniences; neither has a default threshold. Because they are estimated from the current input batch, they are not recommended for a single sample or a very small batch.

```bash
Rscript R/main.R \
  --input /path/to/cohort \
  --genome-build hg19 \
  --output-dir /path/to/results \
  --maf 0.01 \
  --geno 0.05 \
  --plot-mode all
```

Use `--plot-mode none` to disable per-sample plots. Existing outputs are rejected unless `--overwrite` is supplied; overwrite is restricted to the exact output prefix.

## Custom reference PGS

`--reference-pgs FILE` completely replaces the built-in reference distribution. The file must be a tab-delimited table with exactly two columns:

```text
sample_id	pgs_score
sample_001	-4.125
sample_002	-3.718
```

Requirements:

- unique, nonblank `sample_id` values;
- finite numeric `pgs_score` values with nonzero variance;
- fewer than 100 samples: fatal error;
- 100–199 samples: warning;
- absolute adjusted Fisher-Pearson skewness `|G1| > 1`: warning.

LN-PriGI-PRS does not build a reference cohort or calculate its PGS. The supplied reference scores must have been calculated with the same frozen model and scoring convention.

Percentiles use midranks:

```text
100 × (number below + 0.5 × number equal) / reference sample count
```

## Quality control and outputs

Nonzero-weight variant coverage is calculated separately for each sample. Coverage below 90% produces a warning but does not suppress the PGS or percentile. The warning appears in the run report, result table, and plot. Absolute-weight coverage is reported descriptively.

For output prefix `ln_prs`, the program writes:

- `ln_prs.results.tsv`: one row per FID/IID sample;
- `ln_prs.run_report.txt`: configuration, warnings, model/reference summaries, commands, and merged PLINK logs;
- `ln_prs.plots/`: one English PNG per sample when plotting is enabled.

All program-generated reports, tables, plots, messages, and warning codes are in English.

## Frozen assets and reproducibility

The model, built-in reference distribution, fixed 10-sample acceptance dataset, and expected PGS values are versioned with SHA-256 checksums. Runtime assets are verified before scoring. The model is stored as an xz-compressed RDS; fixed PLINK test data are stored in a Deflate ZIP archive.

Run the acceptance test after installing dependencies:

```bash
Rscript tests/run_tests.R
```

If PLINK is not at the test default location, set `LN_PRS_TEST_PLINK` first. When the model changes, keep the fixed test genotypes unchanged, rebuild the expected PGS values, refresh checksums, and rerun the acceptance test. `scripts/build_assets.R` documents the reproducible asset-build inputs.

## License

MIT
