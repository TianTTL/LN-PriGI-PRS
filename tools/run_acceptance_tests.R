#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), winslash = "/", mustWork = TRUE)
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
    rscript <- file.path(
        R.home("bin"),
        if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
}
plink_setting <- Sys.getenv("LN_PRS_TEST_PLINK", unset = "plink")
plink <- if (file.exists(plink_setting)) {
    normalizePath(plink_setting, mustWork = TRUE)
} else {
    Sys.which(plink_setting)
}
if (!nzchar(plink)) {
    stop("PLINK was not found. Set LN_PRS_TEST_PLINK or add plink to PATH.", call. = FALSE)
}

work <- tempfile("ln_prigi_prs_test_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
unzip(file.path(root, "data/test_data/dummy_genome_plink1.zip"), exdir = work)
test_prefix <- file.path(work, "dummy_genome")

output_file <- file.path(work, "results.tsv")
status <- system2(
    rscript,
    c(
        file.path(root, "R/main.R"),
        "--input", test_prefix,
        "--output", output_file,
        "--plot",
        "--plink", plink
    )
)
stopifnot(identical(status, 0L))

observed <- fread(output_file)
expected <- fread(file.path(root, "data/test_data/expected_pgs.tsv"))
merged <- merge(observed, expected, by = c("FID", "IID"))
required_columns <- c(
    "FID", "IID", "PGS", "PERCENTILE", "MODEL_SNPS", "SNPS_USED",
    "VARIANT_COVERAGE", "ABS_WEIGHT_COVERAGE", "REFERENCE_SOURCE"
)
stopifnot(
    nrow(merged) == 10L,
    all(required_columns %in% names(observed)),
    max(abs(merged$PGS - merged$EXPECTED_PGS)) <= 1e-5,
    length(list.files(file.path(work, "results.plots"), pattern = "[.]png$")) == 10L
)

bim_ids <- fread(
    paste0(test_prefix, ".bim"),
    select = 2L,
    header = FALSE
)[[1L]]
model <- as.data.table(readRDS(file.path(root, "data/model/gwas_prior_union.model.rds")))
test_model <- model[MODEL_ID %chin% bim_ids][order(-abs(WEIGHT))]
extract_file <- file.path(work, "ped_smoke.extract")
fwrite(
    test_model[seq_len(1000L), .(SNP_ID = MODEL_ID)],
    extract_file,
    col.names = FALSE
)
ped_prefix <- file.path(work, "ped_smoke")
plink_status <- system2(
    plink,
    c("--bfile", test_prefix, "--extract", extract_file, "--recode", "--out", ped_prefix),
    stdout = FALSE,
    stderr = FALSE
)
stopifnot(identical(plink_status, 0L))
ped <- fread(paste0(ped_prefix, ".ped"), header = FALSE)
missing_genotype_columns <- names(ped)[7:206]
set(
    ped,
    i = 1L,
    j = missing_genotype_columns,
    value = as.list(rep("0", length(missing_genotype_columns)))
)
fwrite(ped, paste0(ped_prefix, ".ped"), sep = " ", col.names = FALSE)
ped_log <- file.path(work, "ped_run.log")
ped_status <- system2(
    rscript,
    c(
        file.path(root, "R/main.R"),
        "--input", ped_prefix,
        "--output", file.path(work, "ped_results.tsv"),
        "--plink", plink
    ),
    stdout = ped_log,
    stderr = ped_log
)
if (!identical(ped_status, 0L)) {
    stop(paste(readLines(ped_log, warn = FALSE), collapse = "\n"), call. = FALSE)
}
ped_results <- fread(file.path(work, "ped_results.tsv"))
first_sample <- ped_results[IID == "test_S1"]
second_sample <- ped_results[IID == "test_S2"]
stopifnot(
    identical(ped_status, 0L),
    nrow(ped_results) == 10L,
    first_sample$VARIANT_COVERAGE < second_sample$VARIANT_COVERAGE,
    first_sample$ABS_WEIGHT_COVERAGE < second_sample$ABS_WEIGHT_COVERAGE
)

small_reference <- file.path(work, "small_reference.tsv")
fwrite(
    data.table(
        sample_id = sprintf("r%03d", 1:99),
        pgs_score = seq_len(99)
    ),
    small_reference,
    sep = "\t"
)
small_reference_status <- system2(
    rscript,
    c(
        file.path(root, "R/main.R"),
        "--input", ped_prefix,
        "--output", file.path(work, "small_reference_results.tsv"),
        "--reference-pgs", small_reference,
        "--plink", plink
    ),
    stdout = FALSE,
    stderr = FALSE
)
stopifnot(!identical(small_reference_status, 0L))

reference_with_extra_column <- file.path(work, "reference_with_extra_column.tsv")
fwrite(
    data.table(
        sample_id = sprintf("r%03d", 1:100),
        pgs_score = seq(-10, 10, length.out = 100),
        population = "test"
    ),
    reference_with_extra_column,
    sep = "\t"
)
extra_column_status <- system2(
    rscript,
    c(
        file.path(root, "R/main.R"),
        "--input", ped_prefix,
        "--output", file.path(work, "extra_column_results.tsv"),
        "--reference-pgs", reference_with_extra_column,
        "--plink", plink
    ),
    stdout = FALSE,
    stderr = FALSE
)
stopifnot(identical(extra_column_status, 0L))

cat("All LN-PriGI-PRS tests passed.\n")
