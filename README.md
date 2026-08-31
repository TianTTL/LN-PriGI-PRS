# LN-PriGI-PRS

LN-PriGI-PRS is a small R command-line tool for calculating a lupus-nephritis polygenic score (PGS) and locating each sample within a reference PGS distribution.

The score and percentile are intended for research. They are not a clinical diagnosis, prognosis, or treatment recommendation.

## Requirements

- R 4.3 or later
- R packages: `data.table` and `dplyr`
- PLINK 1.9 available on `PATH`, or supplied with `--plink`
- Genotypes aligned to hg19/GRCh37

The tool accepts PLINK 1 BED/BIM/FAM or PED/MAP files. It does not accept PLINK 2 PGEN files or VCF directly and does not perform liftover or genotype imputation.

## Usage

Provide the PLINK file prefix without an extension:

```bash
Rscript R/main.R \
  --input /path/to/cohort \
  --output /path/to/results.tsv
```

Specify PLINK explicitly when it is not on `PATH`:

```bash
Rscript R/main.R \
  --input /path/to/cohort \
  --output /path/to/results.tsv \
  --plink /path/to/plink
```

Optional arguments:

```text
--reference-pgs FILE   Use a custom reference PGS table
--plot                Create one reference-distribution plot per sample
--overwrite           Replace the requested output and plot directory
--help                Show command-line help
--version             Show the software version
```

The program does not apply MAF or genotype-missingness filters. Perform any cohort-specific genotype QC before running LN-PriGI-PRS.

## Scoring

Variants are matched using hg19 chromosome, position, and alleles rather than relying on input variant IDs. Non-palindromic strand complements are supported. A/T and C/G palindromic SNPs are not included in the model because their strand direction cannot be resolved reliably across independently prepared genotype datasets.

PGS is the uncentered, unstandardized PLINK `--score ... header sum` result:

```text
PGS = sum(effect-allele dosage × weight)
```

A model variant absent from the input contributes zero. A missing genotype at a matched model variant also contributes zero. Scores are not divided by the number of available variants, so coverage should be considered when comparing datasets.

## Reference distribution and percentile

The built-in reference contains PGS values from 171 control samples scored with the same model and missing-genotype convention used by the tool.

A custom reference must be a tab-delimited file containing at least these columns:

```text
sample_id	pgs_score
sample_001	-7.124
sample_002	-6.835
```

Additional columns are allowed and ignored. Sample IDs must be unique and nonblank, and PGS values must be finite with nonzero variance.

A reference with fewer than 100 samples is rejected. Empirical percentiles based on very small reference sets have coarse resolution and are strongly affected by individual observations, so they are not suitable for the intended comparison.

The tool also calculates adjusted Fisher–Pearson skewness, `G1`. When `|G1| > 1`, it reports a warning because a strongly skewed reference distribution can make percentile differences uneven across the PGS range. The calculation still continues because percentiles are derived from the empirical distribution and do not require normality.

Percentiles use midranks:

```text
100 × (number below + 0.5 × number equal) / reference sample count
```

A percentile describes position within the selected reference distribution; it is not a disease probability.

## Output

The output TSV contains one row per sample:

| Column | Meaning |
| --- | --- |
| `FID`, `IID` | PLINK family and individual identifiers |
| `PGS` | Raw polygenic score |
| `PERCENTILE` | Midrank percentile in the reference distribution |
| `MODEL_SNPS` | Number of variants in the bundled model |
| `SNPS_USED` | Matched, nonmissing model variants for the sample |
| `VARIANT_COVERAGE` | `SNPS_USED / MODEL_SNPS` |
| `ABS_WEIGHT_COVERAGE` | Fraction of total absolute model weight represented by matched, nonmissing variants |
| `REFERENCE_SOURCE` | Built-in or custom reference source |

With `--plot`, plots are written beside the result table in `RESULTS.plots/`. Plotting is disabled by default.

Existing output is not overwritten unless `--overwrite` is supplied.

## Optional installation check

The bundled dummy genotype data can be used to confirm that R, PLINK, and the tool work together:

```bash
Rscript tools/run_acceptance_tests.R
```

If PLINK is not on `PATH`, set `LN_PRS_TEST_PLINK` to its executable path.

## License

MIT
