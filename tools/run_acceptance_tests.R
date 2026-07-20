#!/usr/bin/env Rscript
# Maintainer/developer tool: run the release acceptance suite against the fixed
# 10-sample PLINK cohort. Routine users do not need this script.
#
# Usage:
#   Rscript tools/run_acceptance_tests.R
# PLINK 1.9 is resolved from LN_PRS_TEST_PLINK first, then from PATH. The suite
# verifies asset hashes, golden PGS values, reports, plots, PED/MAP support, and
# strict rejection of a reference distribution with fewer than 100 samples.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), winslash = "/", mustWork = TRUE)
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
plink <- Sys.getenv("LN_PRS_TEST_PLINK", "")
if (!nzchar(plink)) plink <- Sys.which("plink")
if (!nzchar(plink)) stop("PLINK 1.9 was not found. Set LN_PRS_TEST_PLINK or add plink to PATH.", call. = FALSE)
work <- tempfile("ln_prs_test_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
# Verify both the archive and its extracted PLINK files before scoring.
unzip(file.path(root, "data/test_data/dummy_genome_plink1.zip"), exdir = work)
archive_manifest <- fread(file.path(root, "data/test_data/archive_contents.sha256.tsv"))
stopifnot(all(vapply(seq_len(nrow(archive_manifest)), function(i) {
  digest::digest(file.path(work, archive_manifest$path[[i]]), algo = "sha256", file = TRUE) == archive_manifest$sha256[[i]]
}, logical(1))))
out <- file.path(work, "output")
dir.create(out)
cmd <- c(file.path(root, "R/main.R"), "--input", file.path(work, "dummy_genome"),
         "--genome-build", "hg19", "--output-dir", out, "--output-prefix", "acceptance",
         "--plot-mode", "all", "--plink", plink)
status <- system2(rscript, cmd, stdout = file.path(work, "stdout.txt"), stderr = file.path(work, "stderr.txt"))
if (!identical(status, 0L)) stop(paste(readLines(file.path(work, "stderr.txt"), warn = FALSE), collapse = "\n"))
observed <- fread(file.path(out, "acceptance.results.tsv"))
expected <- fread(file.path(root, "data/test_data/expected_pgs.tsv"))
merged <- merge(observed, expected, by = c("FID", "IID"))
stopifnot(nrow(merged) == 10L)
required_columns <- c("FID", "IID", "PGS", "PERCENTILE", "NONZERO_MODEL_SNPS", "NONZERO_SNPS_USED", "NONZERO_VARIANT_COVERAGE", "ABS_WEIGHT_COVERAGE", "QC_STATUS", "PERCENTILE_RELIABLE", "WARNING_CODE", "WARNING_MESSAGE", "REFERENCE_SOURCE")
stopifnot(all(required_columns %in% names(observed)))
stopifnot(max(abs(merged$PGS - merged$EXPECTED_PGS)) <= 1e-5)
stopifnot(all(observed$GENOME_BUILD == "hg19"))
stopifnot(all(observed$NONZERO_VARIANT_COVERAGE >= 0.90))
stopifnot(all(grepl("REFERENCE_SAMPLE_SIZE_LOW", observed$WARNING_CODE)))
stopifnot(file.exists(file.path(out, "acceptance.run_report.txt")))
stopifnot(length(list.files(file.path(out, "acceptance.plots"), pattern = "[.]png$")) == 10L)
report <- readLines(file.path(out, "acceptance.run_report.txt"), warn = FALSE)
stopifnot(any(grepl("REFERENCE_SAMPLE_SIZE_LOW", report, fixed = TRUE)))
# PED/MAP smoke test on a small deterministic variant subset.
bim_ids <- fread(file.path(work, "dummy_genome.bim"), select = 2L, header = FALSE)[[1L]]
extract_file <- file.path(work, "ped_smoke.extract")
fwrite(data.table(SNP_ID = head(bim_ids, 1000L)), extract_file, col.names = FALSE)
ped_prefix <- file.path(work, "ped_smoke")
plink_status <- system2(plink, c("--bfile", file.path(work, "dummy_genome"), "--extract", extract_file, "--recode", "--out", ped_prefix), stdout = FALSE, stderr = FALSE)
stopifnot(identical(plink_status, 0L))
ped_out <- file.path(work, "ped_output"); dir.create(ped_out)
ped_status <- system2(rscript, c(file.path(root, "R/main.R"), "--input", ped_prefix, "--genome-build", "hg19", "--output-dir", ped_out, "--output-prefix", "ped_smoke", "--plot-mode", "none", "--plink", plink), stdout = FALSE, stderr = FALSE)
stopifnot(identical(ped_status, 0L), nrow(fread(file.path(ped_out, "ped_smoke.results.tsv"))) == 10L)

# A custom reference with fewer than 100 samples must terminate strictly.
small_reference <- file.path(work, "small_reference.tsv")
fwrite(data.table(sample_id = sprintf("r%03d", 1:99), pgs_score = seq_len(99)), small_reference, sep = "\t")
small_ref_status <- system2(rscript, c(file.path(root, "R/main.R"), "--input", ped_prefix, "--genome-build", "hg19", "--output-dir", file.path(work, "small_ref_output"), "--reference-pgs", small_reference, "--plot-mode", "none", "--plink", plink), stdout = FALSE, stderr = FALSE)
stopifnot(!identical(small_ref_status, 0L))

cat("All acceptance tests passed.\n")
