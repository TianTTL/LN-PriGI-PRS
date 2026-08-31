#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

VERSION <- "1.0.0"
args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
ROOT <- normalizePath(file.path(dirname(script_arg), ".."), winslash = "/", mustWork = TRUE)

fail <- function(...) stop(paste0(...), call. = FALSE)
"%||%" <- function(x, y) if (is.null(x)) y else x
read_table <- function(path, ...) fread(path, data.table = FALSE, ...)
write_table <- function(x, path, header = TRUE) {
    fwrite(x, path, sep = "\t", quote = FALSE, na = "NA", col.names = header)
}

parse_cli <- function(arguments) {
    flags <- c("--plot", "--overwrite", "--help", "--version")
    valued <- c("--input", "--output", "--reference-pgs", "--plink")
    known <- c(flags, valued)
    parsed <- list()
    i <- 1L
    while (i <= length(arguments)) {
        key <- arguments[[i]]
        if (!key %in% known) fail("Unknown option: ", key)
        name <- substring(key, 3L)
        if (name %in% names(parsed)) fail("Option supplied more than once: ", key)
        if (key %in% flags) {
            parsed[[name]] <- TRUE
            i <- i + 1L
        } else {
            if (i == length(arguments) || startsWith(arguments[[i + 1L]], "--")) {
                fail("Missing value for ", key)
            }
            parsed[[name]] <- arguments[[i + 1L]]
            i <- i + 2L
        }
    }
    parsed
}

usage <- function() {
    cat(paste0(
        "LN-PriGI-PRS ", VERSION, "\n\n",
        "Usage:\n  Rscript R/main.R --input PREFIX --output RESULTS.tsv [options]\n\n",
        "Required:\n",
        "  --input PREFIX           PLINK BED/BIM/FAM or PED/MAP prefix\n",
        "  --output FILE            Output TSV file\n\n",
        "Options:\n",
        "  --reference-pgs FILE     Reference TSV containing sample_id and pgs_score\n",
        "  --plink PATH             PLINK 1.9 executable [search PATH]\n",
        "  --plot                    Create one plot per sample\n",
        "  --overwrite               Replace the requested output and plot directory\n",
        "  --help                    Show help\n",
        "  --version                 Show version\n"
    ))
}

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
    usage()
    quit(status = 0L)
}
if (isTRUE(opt$version)) {
    cat(VERSION, "\n")
    quit(status = 0L)
}
missing_opt <- setdiff(c("input", "output"), names(opt))
if (length(missing_opt)) {
    fail("Missing required option(s): --", paste(missing_opt, collapse = ", --"))
}

input_prefix <- normalizePath(opt$input, winslash = "/", mustWork = FALSE)
has_file <- function(extension) file.exists(paste0(input_prefix, extension))
has_binary_set <- all(vapply(c(".bed", ".bim", ".fam"), has_file, logical(1)))
has_text_set <- all(vapply(c(".ped", ".map"), has_file, logical(1)))
any_binary_file <- any(vapply(c(".bed", ".bim", ".fam"), has_file, logical(1)))
any_text_file <- any(vapply(c(".ped", ".map"), has_file, logical(1)))
if (has_binary_set && has_text_set) fail("Both BED/BIM/FAM and PED/MAP were found.")
if (!has_binary_set && !has_text_set) {
    fail(if (any_binary_file || any_text_file) {
        "The PLINK input file set is incomplete."
    } else {
        "No supported PLINK input files were found."
    })
}

plink_setting <- opt$plink %||% "plink"
plink <- if (file.exists(plink_setting)) {
    normalizePath(plink_setting, winslash = "/", mustWork = TRUE)
} else {
    Sys.which(plink_setting)
}
if (!nzchar(plink)) {
    fail("PLINK was not found. Add PLINK 1.9 to PATH or use --plink PATH.")
}

output_path <- normalizePath(opt$output, winslash = "/", mustWork = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dirname(output_path))) fail("Cannot create output directory.")
plot_dir <- file.path(
    dirname(output_path),
    paste0(tools::file_path_sans_ext(basename(output_path)), ".plots")
)
existing <- output_path[file.exists(output_path)]
if (isTRUE(opt$plot) && dir.exists(plot_dir)) existing <- c(existing, plot_dir)
if (length(existing) && !isTRUE(opt$overwrite)) {
    fail("Output already exists. Use --overwrite to replace: ", existing[[1]])
}
if (isTRUE(opt$overwrite)) {
    if (file.exists(output_path)) unlink(output_path, force = TRUE)
    if (isTRUE(opt$plot) && dir.exists(plot_dir)) {
        unlink(plot_dir, recursive = TRUE, force = TRUE)
    }
}

work <- tempfile("ln_prigi_prs_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

run_plink <- function(label, arguments) {
    log_path <- file.path(work, paste0(label, ".log"))
    status <- system2(plink, arguments, stdout = log_path, stderr = log_path)
    log_text <- readLines(log_path, warn = FALSE)
    if (!identical(status, 0L)) {
        fail("PLINK failed during ", label, ".\n", paste(tail(log_text, 20L), collapse = "\n"))
    }
    invisible(log_text)
}

model <- readRDS(file.path(ROOT, "data/model/gwas_prior_union.model.rds")) |>
    as_tibble() |>
    mutate(
        EFFECT_ALLELE = toupper(.data$EFFECT_ALLELE),
        OTHER_ALLELE = toupper(.data$OTHER_ALLELE)
    )
model_n <- nrow(model)
model_abs_weight <- sum(abs(model$WEIGHT))

base_prefix <- input_prefix
if (has_text_set) {
    base_prefix <- file.path(work, "ped_converted")
    run_plink("convert_ped", c("--file", input_prefix, "--make-bed", "--out", base_prefix))
}

bim <- read_table(
    paste0(base_prefix, ".bim"),
    col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2")
) |>
    mutate(
        CHR = suppressWarnings(as.integer(.data$CHR)),
        POS = suppressWarnings(as.integer(.data$POS)),
        A1 = toupper(.data$A1),
        A2 = toupper(.data$A2)
    )
duplicate_coordinates <- bim |>
    count(.data$CHR, .data$POS, name = "n") |>
    filter(.data$n > 1L) |>
    select("CHR", "POS")
if (nrow(duplicate_coordinates)) {
    bim <- anti_join(bim, duplicate_coordinates, by = c("CHR", "POS"))
}
duplicate_ids <- bim |>
    count(.data$SNP_ID, name = "n") |>
    filter(.data$n > 1L) |>
    select("SNP_ID")
if (nrow(duplicate_ids)) bim <- anti_join(bim, duplicate_ids, by = "SNP_ID")

complement <- function(x) chartr("ACGT", "TGCA", x)
matched <- bim |>
    transmute(
        SNP_ID = .data$SNP_ID,
        CHR = .data$CHR,
        POS = .data$POS,
        USER_A1 = .data$A1,
        USER_A2 = .data$A2
    ) |>
    inner_join(model, by = c("CHR", "POS")) |>
    mutate(
        DIRECT_MATCH = (
            .data$EFFECT_ALLELE == .data$USER_A1 &
                .data$OTHER_ALLELE == .data$USER_A2
        ) | (
            .data$EFFECT_ALLELE == .data$USER_A2 &
                .data$OTHER_ALLELE == .data$USER_A1
        ),
        COMPLEMENT_MATCH = (
            complement(.data$EFFECT_ALLELE) == .data$USER_A1 &
                complement(.data$OTHER_ALLELE) == .data$USER_A2
        ) | (
            complement(.data$EFFECT_ALLELE) == .data$USER_A2 &
                complement(.data$OTHER_ALLELE) == .data$USER_A1
        )
    ) |>
    filter(.data$DIRECT_MATCH | .data$COMPLEMENT_MATCH) |>
    mutate(
        USER_EFFECT = if_else(.data$DIRECT_MATCH, .data$EFFECT_ALLELE, complement(.data$EFFECT_ALLELE)),
        USER_OTHER = if_else(.data$USER_EFFECT == .data$USER_A1, .data$USER_A2, .data$USER_A1)
    )
if (!nrow(matched)) fail("No model variants matched the input genotype data.")

score_file <- file.path(work, "score.tsv")
a2_file <- file.path(work, "a2.tsv")
extract_file <- file.path(work, "extract.tsv")
write_table(select(matched, "SNP_ID", "USER_EFFECT", "WEIGHT"), score_file)
write_table(select(matched, "SNP_ID", "USER_OTHER"), a2_file, header = FALSE)
write_table(select(matched, "SNP_ID"), extract_file, header = FALSE)

selected_prefix <- file.path(work, "selected")
run_plink(
    "extract",
    c("--bfile", base_prefix, "--extract", extract_file, "--make-bed", "--out", selected_prefix)
)
selected_bim <- read_table(
    paste0(selected_prefix, ".bim"),
    col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2")
)
retained <- semi_join(matched, select(selected_bim, "SNP_ID"), by = "SNP_ID")
if (!nrow(retained)) fail("No model variants remained after PLINK extraction.")

# Measure missingness before recoding so coverage reflects observed genotypes;
# missing calls are then set to the non-effect allele so they contribute zero to PGS.
missing_prefix <- file.path(work, "missing")
run_plink("missingness", c("--bfile", selected_prefix, "--missing", "--out", missing_prefix))
missingness <- read_table(paste0(missing_prefix, ".imiss"))
if (anyDuplicated(select(missingness, "FID", "IID"))) fail("FID/IID combinations must be unique.")
retained_n <- nrow(retained)
matched_abs_weight <- sum(abs(retained$WEIGHT))
missingness <- missingness |>
    mutate(
        SNPS_USED = pmax(0, retained_n - .data$N_MISS),
        VARIANT_COVERAGE = .data$SNPS_USED / model_n
    )

missing_variant_ids <- read_table(paste0(missing_prefix, ".lmiss")) |>
    filter(.data$N_MISS > 0L) |>
    pull(.data$SNP)
if (length(missing_variant_ids)) {
    missing_variant_file <- file.path(work, "missing_variants.tsv")
    write_table(data.frame(SNP_ID = missing_variant_ids), missing_variant_file, header = FALSE)
    missing_genotype_prefix <- file.path(work, "missing_genotype_data")
    run_plink(
        "missing_genotypes",
        c(
            "--bfile", selected_prefix,
            "--extract", missing_variant_file,
            "--recode", "A",
            "--out", missing_genotype_prefix
        )
    )
    missing_genotypes <- read_table(paste0(missing_genotype_prefix, ".raw"))
    genotype_columns <- names(missing_genotypes)[-(1:6)]
    genotype_variant_ids <- sub("_[^_]+$", "", genotype_columns)
    variant_abs_weights <- setNames(abs(retained$WEIGHT), retained$SNP_ID)
    column_abs_weights <- unname(variant_abs_weights[genotype_variant_ids])
    if (anyNA(column_abs_weights)) fail("Could not match missing-genotype columns to model weights.")
    missing_abs_weight <- as.numeric(
        is.na(as.matrix(missing_genotypes[, genotype_columns, drop = FALSE])) %*%
            column_abs_weights
    )
    sample_abs_weight_coverage <- missing_genotypes |>
        transmute(
            FID = .data$FID,
            IID = .data$IID,
            ABS_WEIGHT_COVERAGE = (matched_abs_weight - missing_abs_weight) / model_abs_weight
        )
    missingness <- missingness |>
        inner_join(sample_abs_weight_coverage, by = c("FID", "IID"))
} else {
    missingness <- missingness |>
        mutate(ABS_WEIGHT_COVERAGE = matched_abs_weight / model_abs_weight)
}

score_ready_prefix <- file.path(work, "score_ready")
a2_log <- run_plink(
    "encode_missing",
    c(
        "--bfile", selected_prefix,
        "--a2-allele", a2_file, "2", "1",
        "--fill-missing-a2",
        "--make-bed",
        "--out", score_ready_prefix
    )
)
if (any(grepl("Impossible A2 allele assignment", a2_log, fixed = TRUE))) {
    fail("Allele alignment failed during PLINK A2 assignment.")
}

pgs_prefix <- file.path(work, "pgs")
run_plink(
    "score",
    c(
        "--bfile", score_ready_prefix,
        "--score", score_file, "1", "2", "3", "header", "sum",
        "--out", pgs_prefix
    )
)
pgs <- read_table(paste0(pgs_prefix, ".profile")) |>
    transmute(FID = .data$FID, IID = .data$IID, PGS = .data$SCORESUM)

adjusted_skewness <- function(x) {
    n <- length(x)
    n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / sd(x))^3)
}
read_reference <- function(path = NULL) {
    source_path <- if (is.null(path)) {
        file.path(ROOT, "data/reference/reference_pgs.tsv")
    } else {
        normalizePath(path, winslash = "/", mustWork = TRUE)
    }
    reference <- read_table(source_path, sep = "\t")
    missing_columns <- setdiff(c("sample_id", "pgs_score"), names(reference))
    if (length(missing_columns)) {
        fail("Reference PGS is missing column(s): ", paste(missing_columns, collapse = ", "))
    }
    reference <- reference |>
        transmute(
            sample_id = trimws(as.character(.data$sample_id)),
            pgs_score = suppressWarnings(as.numeric(.data$pgs_score))
        )
    if (
        anyNA(reference$sample_id) ||
            any(!nzchar(reference$sample_id)) ||
            anyDuplicated(reference$sample_id) ||
            any(!is.finite(reference$pgs_score))
    ) {
        fail("Reference sample IDs must be unique and nonblank, and PGS values must be finite.")
    }
    if (nrow(reference) < 100L) fail("Reference PGS must contain at least 100 samples.")
    if (sd(reference$pgs_score) == 0) fail("Reference PGS must have nonzero variance.")

    skewness <- adjusted_skewness(reference$pgs_score)
    if (abs(skewness) > 1) {
        warning(
            "Reference PGS is skewed (adjusted Fisher-Pearson G1 = ",
            round(skewness, 3),
            "). Percentiles are still calculated from the empirical distribution.",
            call. = FALSE
        )
    }
    list(
        data = reference,
        source = if (is.null(path)) "built-in" else source_path,
        skewness = skewness
    )
}

reference_info <- read_reference(opt[["reference-pgs"]])
reference <- reference_info$data
midrank_percentile <- function(score, reference_scores) {
    100 * (sum(reference_scores < score) + 0.5 * sum(reference_scores == score)) /
        length(reference_scores)
}

results <- pgs |>
    inner_join(
        select(missingness, "FID", "IID", "SNPS_USED", "VARIANT_COVERAGE", "ABS_WEIGHT_COVERAGE"),
        by = c("FID", "IID")
    ) |>
    mutate(
        FID = if_else(as.character(.data$FID) == "0", NA_character_, as.character(.data$FID)),
        IID = as.character(.data$IID),
        PERCENTILE = vapply(
            .data$PGS,
            midrank_percentile,
            numeric(1),
            reference_scores = reference$pgs_score
        ),
        MODEL_SNPS = model_n,
        REFERENCE_SOURCE = reference_info$source
    ) |>
    transmute(
        FID = .data$FID,
        IID = .data$IID,
        PGS = round(.data$PGS, 6),
        PERCENTILE = round(.data$PERCENTILE, 2),
        MODEL_SNPS = .data$MODEL_SNPS,
        SNPS_USED = .data$SNPS_USED,
        VARIANT_COVERAGE = round(.data$VARIANT_COVERAGE, 4),
        ABS_WEIGHT_COVERAGE = round(.data$ABS_WEIGHT_COVERAGE, 4),
        REFERENCE_SOURCE = .data$REFERENCE_SOURCE
    )
write_table(results, output_path)

if (isTRUE(opt$plot)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- function(x) {
        x <- gsub("[^A-Za-z0-9._-]+", "_", x)
        ifelse(nzchar(x), x, "sample")
    }
    for (i in seq_len(nrow(results))) {
        sample_label <- paste(na.omit(c(results$FID[[i]], results$IID[[i]])), collapse = "__")
        plot_path <- file.path(plot_dir, sprintf("%03d_%s.png", i, safe_name(sample_label)))
        png(plot_path, width = 1200, height = 800, res = 140)
        hist(
            reference$pgs_score,
            breaks = "FD",
            col = "#D9E6F2",
            border = "white",
            main = paste0("LN-PriGI-PRS: ", sample_label),
            xlab = "PGS"
        )
        abline(v = results$PGS[[i]], col = "#C43C39", lwd = 3)
        mtext(
            paste0(
                "PGS = ", results$PGS[[i]],
                " | percentile = ", results$PERCENTILE[[i]],
                " | SNP coverage = ", sprintf("%.1f%%", 100 * results$VARIANT_COVERAGE[[i]])
            ),
            side = 3,
            line = 0.2
        )
        dev.off()
    }
}

message(
    "LN-PriGI-PRS completed.\n",
    "  Samples: ", nrow(results), "\n",
    "  Model variants: ", model_n, "\n",
    "  Matched variants: ", retained_n, "\n",
    "  Reference samples: ", nrow(reference), "\n",
    "  Reference skewness G1: ", round(reference_info$skewness, 3), "\n",
    "  Results: ", output_path,
    if (isTRUE(opt$plot)) paste0("\n  Plots: ", plot_dir) else ""
)
