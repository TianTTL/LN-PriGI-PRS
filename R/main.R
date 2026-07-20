#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(digest)
    library(dplyr)
})

# Resolve the repository root so bundled assets can be found from any working directory.
VERSION <- "0.1.0"
args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
ROOT <- normalizePath(file.path(dirname(script_arg), ".."), winslash = "/", mustWork = TRUE)

fail <- function(...) {
    stop(paste0(...), call. = FALSE)
}

"%||%" <- function(x, y) {
    if (is.null(x)) y else x
}

read_table <- function(path, ...) {
    data.table::fread(path, data.table = FALSE, ...)
}

write_table <- function(data, path, column_names = TRUE) {
    data.table::fwrite(
        data,
        path,
        sep = "\t",
        quote = FALSE,
        na = "NA",
        col.names = column_names
    )
}

# Parse the deliberately small CLI without adding another runtime dependency.
parse_cli <- function(arguments) {
    flags <- c("--overwrite", "--help", "--version")
    parsed <- list()
    index <- 1L

    while (index <= length(arguments)) {
        key <- arguments[[index]]
        if (!startsWith(key, "--")) {
            fail("Unexpected argument: ", key)
        }
        if (key %in% flags) {
            parsed[[substring(key, 3L)]] <- TRUE
            index <- index + 1L
            next
        }
        if (index == length(arguments) || startsWith(arguments[[index + 1L]], "--")) {
            fail("Missing value for ", key)
        }
        parsed[[substring(key, 3L)]] <- arguments[[index + 1L]]
        index <- index + 2L
    }

    parsed
}

usage <- function() {
    cat(paste0(
        "LN-PriGI-PRS ", VERSION, "\n\n",
        "Usage:\n  Rscript R/main.R --input PREFIX --genome-build hg19 --output-dir DIR [options]\n\n",
        "Required:\n  --input PATH              PLINK 1 BED/BIM/FAM or PED/MAP prefix\n",
        "  --genome-build hg19      Explicit genome-build declaration; only hg19 is accepted\n",
        "  --output-dir DIR         Output directory\n\n",
        "Options:\n  --output-prefix NAME      Output prefix [ln_prs]\n",
        "  --reference-pgs FILE     Replacement reference TSV: sample_id, pgs_score\n",
        "  --plink PATH             PLINK 1.9 executable [search PATH]\n",
        "  --maf FLOAT              Optional PLINK MAF filter; no default\n",
        "  --geno FLOAT             Optional PLINK missingness filter; no default\n",
        "  --plot-mode all|none     Per-sample plots [all]\n",
        "  --overwrite              Replace this exact output prefix\n",
        "  --help                    Show help\n",
        "  --version                 Show version\n"
    ))
}

# Validate the public command-line contract before reading large genotype files.
opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
    usage()
    quit(status = 0L)
}
if (isTRUE(opt$version)) {
    cat(VERSION, "\n")
    quit(status = 0L)
}
if ("input-prefix" %in% names(opt) && !"input" %in% names(opt)) {
    opt$input <- opt[["input-prefix"]]
}

required <- c("input", "genome-build", "output-dir")
missing_opt <- required[!required %in% names(opt)]
if (length(missing_opt)) {
    fail("Missing required option(s): --", paste(missing_opt, collapse = ", --"))
}
if (!identical(opt[["genome-build"]], "hg19")) {
    fail("Only hg19 input is accepted. Liftover is not provided.")
}

plot_mode <- opt[["plot-mode"]] %||% "all"
if (!plot_mode %in% c("all", "none")) {
    fail("--plot-mode must be 'all' or 'none'.")
}
out_prefix <- opt[["output-prefix"]] %||% "ln_prs"
if (!grepl("^[A-Za-z0-9._-]+$", out_prefix)) {
    fail("--output-prefix may contain only letters, digits, dot, underscore, and hyphen.")
}

validate_num <- function(name, lower, upper) {
    if (!name %in% names(opt)) {
        return(NULL)
    }
    value <- suppressWarnings(as.numeric(opt[[name]]))
    if (!is.finite(value) || value < lower || value > upper) {
        fail("--", name, " must be between ", lower, " and ", upper, ".")
    }
    value
}

maf <- validate_num("maf", 0, 0.5)
geno <- validate_num("geno", 0, 1)
# Detect exactly one supported PLINK 1 file set and reject partial or ambiguous inputs.
input_prefix <- normalizePath(opt[["input"]], winslash = "/", mustWork = FALSE)
has_input_file <- function(extension) {
    file.exists(paste0(input_prefix, extension))
}

binary <- all(vapply(c(".bed", ".bim", ".fam"), has_input_file, logical(1)))
text <- all(vapply(c(".ped", ".map"), has_input_file, logical(1)))
any_binary <- any(vapply(c(".bed", ".bim", ".fam"), has_input_file, logical(1)))
any_text <- any(vapply(c(".ped", ".map"), has_input_file, logical(1)))
if (binary && text) {
    fail("Both BED/BIM/FAM and PED/MAP were found. Provide an unambiguous prefix.")
}
if (!binary && !text) {
    message <- if (any_binary || any_text) {
        "The PLINK input file set is incomplete."
    } else {
        "No supported PLINK 1 input files were found."
    }
    fail(message)
}
input_mode <- if (binary) "BED/BIM/FAM" else "PED/MAP"

# Resolve and verify PLINK before creating any analysis output.
plink <- opt$plink %||% Sys.which("plink")
if (!nzchar(plink) || !file.exists(plink)) {
    fail("PLINK was not found. Install PLINK 1.9 or provide --plink PATH.")
}
plink <- normalizePath(plink, winslash = "/", mustWork = TRUE)
version_lines <- suppressWarnings(system2(plink, "--version", stdout = TRUE, stderr = TRUE))
if (!any(grepl("PLINK v1\\.9", version_lines))) {
    fail("The selected executable is not PLINK 1.9: ", paste(version_lines, collapse = " "))
}

# Reserve the exact output prefix and isolate all intermediate files in a temporary directory.
out_dir <- normalizePath(opt[["output-dir"]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
    fail("Cannot create output directory: ", out_dir)
}
results_path <- file.path(out_dir, paste0(out_prefix, ".results.tsv"))
report_path <- file.path(out_dir, paste0(out_prefix, ".run_report.txt"))
plots_dir <- file.path(out_dir, paste0(out_prefix, ".plots"))
targets <- c(results_path, report_path, plots_dir)
existing <- targets[file.exists(targets) | dir.exists(targets)]
if (length(existing) && !isTRUE(opt$overwrite)) {
    fail("Output already exists. Use --overwrite to replace this exact prefix: ", existing[[1]])
}
if (isTRUE(opt$overwrite)) {
    for (path in existing) {
        if (dir.exists(path)) {
            unlink(path, recursive = TRUE, force = TRUE)
        } else {
            unlink(path, force = TRUE)
        }
    }
}

work <- tempfile(paste0("ln_prs_", out_prefix, "_"), tmpdir = out_dir)
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
logs <- character()
commands <- character()

# Execute PLINK consistently while collecting commands, exit codes, stdout, and stderr.
run_plink <- function(label, arguments) {
    log_path <- file.path(work, paste0(label, ".log"))
    command_text <- paste(shQuote(plink), paste(shQuote(arguments), collapse = " "))
    status <- system2(plink, arguments, stdout = log_path, stderr = log_path)
    commands <<- c(commands, paste0(command_text, " [exit=", status, "]"))
    log_text <- readLines(log_path, warn = FALSE)
    logs <<- c(logs, paste0("\n===== PLINK: ", label, " =====\n"), log_text)

    if (!identical(status, 0L)) {
        failure_report <- c(
            "LN-PriGI-PRS failed.",
            paste0("Stage: ", label),
            "",
            commands,
            logs
        )
        writeLines(failure_report, report_path, useBytes = TRUE)
        fail("PLINK failed during ", label, ". Details were written to ", report_path)
    }
    log_text
}

# Verify only the frozen assets needed by this run before loading the model.
manifest_path <- file.path(ROOT, "data/manifest.sha256.tsv")
manifest <- read_table(manifest_path)
required_assets <- c(
    "model/gwas_prior_union.model.rds",
    "model/model_metadata.tsv"
)
if (!"reference-pgs" %in% names(opt)) {
    required_assets <- c(
        required_assets,
        "reference/reference_pgs.tsv",
        "reference/reference_metadata.tsv"
    )
}
manifest <- manifest |>
    filter(.data$path %in% required_assets)

for (index in seq_len(nrow(manifest))) {
    asset <- file.path(ROOT, "data", manifest$path[[index]])
    checksum_is_valid <- file.exists(asset) &&
        digest(asset, algo = "sha256", file = TRUE) == manifest$sha256[[index]]
    if (!checksum_is_valid) {
        fail("Asset checksum verification failed: ", manifest$path[[index]])
    }
}

# Load the compact frozen model; data.table is used only for high-volume file I/O.
model_metadata <- read_table(file.path(ROOT, "data/model/model_metadata.tsv"))
model <- readRDS(file.path(ROOT, "data/model/gwas_prior_union.model.rds")) |>
    as_tibble() |>
    mutate(
        EFFECT_ALLELE = toupper(.data$EFFECT_ALLELE),
        OTHER_ALLELE = toupper(.data$OTHER_ALLELE)
    )
model_n <- nrow(model)
model_abs <- sum(abs(model$WEIGHT))

# Convert PED/MAP once so the remaining workflow operates on one PLINK representation.
base_prefix <- input_prefix
if (text) {
    converted <- file.path(work, "ped_converted")
    run_plink("convert_ped", c("--file", input_prefix, "--make-bed", "--out", converted))
    base_prefix <- converted
}

# Remove ambiguous input coordinates and duplicate IDs before model alignment.
bim <- read_table(
    paste0(base_prefix, ".bim"),
    col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2")
) |>
    mutate(
        A1 = toupper(.data$A1),
        A2 = toupper(.data$A2),
        POS = as.integer(.data$POS)
    )
dup_coord <- bim |>
    count(.data$CHR, .data$POS, name = "count") |>
    filter(.data$count > 1L) |>
    select("CHR", "POS")
if (nrow(dup_coord)) {
    bim <- anti_join(bim, dup_coord, by = c("CHR", "POS"))
}
dup_id <- bim |>
    count(.data$SNP_ID, name = "count") |>
    filter(.data$count > 1L) |>
    select("SNP_ID")
if (nrow(dup_id)) {
    bim <- anti_join(bim, dup_id, by = "SNP_ID")
}

# Align variants by hg19 position and the complete biallelic allele set on either strand.
complement_allele <- function(allele) {
    chartr("ACGT", "TGCA", allele)
}

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
            complement_allele(.data$EFFECT_ALLELE) == .data$USER_A1 &
                complement_allele(.data$OTHER_ALLELE) == .data$USER_A2
        ) | (
            complement_allele(.data$EFFECT_ALLELE) == .data$USER_A2 &
                complement_allele(.data$OTHER_ALLELE) == .data$USER_A1
        )
    ) |>
    filter(.data$DIRECT_MATCH | .data$COMPLEMENT_MATCH) |>
    mutate(
        USER_EFFECT = if_else(
            .data$DIRECT_MATCH,
            .data$EFFECT_ALLELE,
            complement_allele(.data$EFFECT_ALLELE)
        ),
        USER_OTHER = if_else(
            .data$USER_EFFECT == .data$USER_A1,
            .data$USER_A2,
            .data$USER_A1
        )
    )
if (!nrow(matched)) {
    fail("No eligible model variants matched the input genotype data.")
}
direct_match_n <- sum(matched$DIRECT_MATCH)
complement_match_n <- sum(matched$COMPLEMENT_MATCH & !matched$DIRECT_MATCH)

# Apply optional cohort-level filters only after allele-compatible model variants are known.
score_initial <- file.path(work, "score_initial.tsv")
matched |>
    select("SNP_ID", "USER_EFFECT", "WEIGHT") |>
    write_table(score_initial)

qc <- file.path(work, "qc")
qc_args <- c("--bfile", base_prefix, "--extract", score_initial)
if (!is.null(maf)) {
    qc_args <- c(qc_args, "--maf", format(maf, scientific = FALSE))
}
if (!is.null(geno)) {
    qc_args <- c(qc_args, "--geno", format(geno, scientific = FALSE))
}
run_plink("extract_and_filter", c(qc_args, "--make-bed", "--out", qc))

# Read the post-filter BIM to retain only variants PLINK kept.
qc_bim <- read_table(
    paste0(qc, ".bim"),
    col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2")
)
retained <- matched |>
    semi_join(select(qc_bim, "SNP_ID"), by = "SNP_ID")
if (!nrow(retained)) {
    fail("No model variants remained after optional filters.")
}

score <- file.path(work, "score.tsv")
a2 <- file.path(work, "a2.tsv")
retained |>
    select("SNP_ID", "USER_EFFECT", "WEIGHT") |>
    write_table(score)
retained |>
    select("SNP_ID", "USER_OTHER") |>
    write_table(a2, column_names = FALSE)

# Measure per-sample missingness before creating the zero-contribution scoring representation.
missing_prefix <- file.path(work, "missing")
run_plink("missingness", c("--bfile", qc, "--missing", "--out", missing_prefix))
imiss <- read_table(paste0(missing_prefix, ".imiss"))
if (anyDuplicated(select(imiss, "FID", "IID"))) {
    fail("FID/IID combinations must be unique.")
}
retained_n <- nrow(retained)
imiss <- imiss |>
    mutate(
        NONZERO_SNPS_USED = pmax(0, retained_n - .data$N_MISS),
        NONZERO_VARIANT_COVERAGE = .data$NONZERO_SNPS_USED / model_n
    )
abs_coverage <- sum(abs(retained$WEIGHT)) / model_abs

# Encode missing target genotypes as the non-effect allele in a temporary copy, yielding zero contribution.
score_ready <- file.path(work, "score_ready")
a2_log <- run_plink(
    "encode_missing_as_zero",
    c(
        "--bfile", qc,
        "--a2-allele", a2, "2", "1",
        "--fill-missing-a2",
        "--make-bed",
        "--out", score_ready
    )
)
if (any(grepl("Impossible A2 allele assignment", a2_log, fixed = TRUE))) {
    fail("Internal allele alignment verification failed during A2 assignment.")
}

prs_prefix <- file.path(work, "prs")
run_plink(
    "score",
    c(
        "--bfile", score_ready,
        "--score", score, "1", "2", "3", "header", "sum",
        "--out", prs_prefix
    )
)
prs <- read_table(paste0(prs_prefix, ".profile")) |>
    transmute(FID = .data$FID, IID = .data$IID, PGS = .data$SCORESUM)
res <- prs |>
    inner_join(
        select(
            imiss,
            "FID",
            "IID",
            "NONZERO_SNPS_USED",
            "NONZERO_VARIANT_COVERAGE"
        ),
        by = c("FID", "IID")
    ) |>
    mutate(FID = if_else(.data$FID == "0", NA_character_, as.character(.data$FID)))

# Validate either the built-in or replacement reference distribution with strict stop conditions.
read_reference <- function(path = NULL) {
    source_path <- if (is.null(path)) {
        file.path(ROOT, "data/reference/reference_pgs.tsv")
    } else {
        path
    }
    reference <- read_table(source_path, sep = "\t")
    if (!identical(names(reference), c("sample_id", "pgs_score"))) {
        fail("Reference PGS must be a tab-delimited table with exactly: sample_id, pgs_score.")
    }
    if (nrow(reference) < 100L) {
        fail(
            "The reference PGS has too few samples to represent population distribution characteristics ",
            "(n < 100)."
        )
    }
    invalid_ids <- anyNA(reference$sample_id) ||
        any(!nzchar(trimws(reference$sample_id))) ||
        anyDuplicated(reference$sample_id)
    if (invalid_ids) {
        fail("Reference sample_id values must be nonblank and unique.")
    }
    reference <- reference |>
        mutate(pgs_score = suppressWarnings(as.numeric(.data$pgs_score)))
    if (any(!is.finite(reference$pgs_score)) || sd(reference$pgs_score) == 0) {
        fail("Reference pgs_score values must be finite and have nonzero variance.")
    }

    sample_count <- nrow(reference)
    standardized <- (reference$pgs_score - mean(reference$pgs_score)) / sd(reference$pgs_score)
    skewness <- sample_count / ((sample_count - 1) * (sample_count - 2)) * sum(standardized^3)
    source_name <- if (is.null(path)) {
        "built-in:LN_WGS_controls"
    } else {
        normalizePath(path, winslash = "/", mustWork = TRUE)
    }
    list(data = reference, source = source_name, skew = skewness)
}

# Calculate empirical midrank percentiles while retaining raw, unstandardized PGS values.
reference <- read_reference(opt[["reference-pgs"]])
reference_values <- reference$data$pgs_score
reference_n <- length(reference_values)
calculate_percentile <- function(value) {
    below <- sum(reference_values < value)
    equal <- sum(reference_values == value)
    round(100 * (below + 0.5 * equal) / reference_n, 2)
}
res <- res |>
    mutate(
        PERCENTILE = vapply(.data$PGS, calculate_percentile, numeric(1)),
        NONZERO_MODEL_SNPS = model_n,
        ABS_WEIGHT_COVERAGE = abs_coverage,
        warnings__ = ""
    )

# Accumulate machine-readable warnings without suppressing sample scores or percentiles.
add_warning <- function(data, code, condition = rep(TRUE, nrow(data))) {
    data |>
        mutate(
            warnings__ = if_else(
                condition,
                if_else(
                    nzchar(.data$warnings__),
                    paste0(.data$warnings__, ";", code),
                    code
                ),
                .data$warnings__
            )
        )
}

if (reference_n < 200L) {
    res <- add_warning(res, "REFERENCE_SAMPLE_SIZE_LOW")
}
if (abs(reference$skew) > 1) {
    res <- add_warning(res, "REFERENCE_DISTRIBUTION_SKEWED")
}
res <- add_warning(
    res,
    "NONZERO_VARIANT_COVERAGE_BELOW_90_PERCENT",
    res$NONZERO_VARIANT_COVERAGE < 0.90
)

# Translate warning codes into stable English messages for tables and reports.
warning_text <- c(
    REFERENCE_SAMPLE_SIZE_LOW = paste(
        "The reference sample size is below 200 and may not fully represent",
        "the population distribution."
    ),
    REFERENCE_DISTRIBUTION_SKEWED = paste(
        "The reference PGS distribution has absolute adjusted Fisher-Pearson",
        "skewness above 1."
    ),
    NONZERO_VARIANT_COVERAGE_BELOW_90_PERCENT = paste(
        "Nonzero-weight model variant coverage is below 90%."
    ),
    PLOT_FAILED = "The per-sample plot could not be generated."
)
warning_message <- function(codes) {
    split_codes <- strsplit(codes, ";", fixed = TRUE)[[1L]]
    paste(unname(warning_text[split_codes[nzchar(split_codes)]]), collapse = " ")
}
res <- res |>
    mutate(
        GENOME_BUILD = "hg19",
        REFERENCE_SOURCE = reference$source,
        PERCENTILE_RELIABLE = .data$NONZERO_VARIANT_COVERAGE >= 0.90,
        QC_STATUS = if_else(nzchar(.data$warnings__), "WARN", "PASS"),
        WARNING_CODE = .data$warnings__,
        WARNING_MESSAGE = vapply(.data$warnings__, warning_message, character(1))
    ) |>
    select(
        all_of(c(
            "FID",
            "IID",
            "PGS",
            "PERCENTILE",
            "NONZERO_MODEL_SNPS",
            "NONZERO_SNPS_USED",
            "NONZERO_VARIANT_COVERAGE",
            "ABS_WEIGHT_COVERAGE",
            "QC_STATUS",
            "PERCENTILE_RELIABLE",
            "WARNING_CODE",
            "WARNING_MESSAGE",
            "REFERENCE_SOURCE",
            "GENOME_BUILD"
        ))
    )

# Finalize tabular output atomically so interrupted runs do not leave a partial results file.
write_atomic <- function(data, path) {
    temporary_path <- paste0(path, ".tmp")
    write_table(data, temporary_path)
    if (!file.rename(temporary_path, path)) {
        fail("Cannot finalize output: ", path)
    }
}

# Render one collision-safe plot per sample and downgrade isolated plot failures to warnings.
if (plot_mode == "all") {
    dir.create(plots_dir)
    sanitize_filename <- function(value) {
        gsub("[^A-Za-z0-9._-]", "_", value)
    }
    plot_stems <- paste0(sanitize_filename(res$FID), "__", sanitize_filename(res$IID))
    collisions <- duplicated(plot_stems) | duplicated(plot_stems, fromLast = TRUE)
    collision_ids <- paste(res$FID[collisions], res$IID[collisions], sep = "\t")
    collision_hashes <- substr(
        vapply(collision_ids, digest, character(1), algo = "sha256", serialize = FALSE),
        1L,
        8L
    )
    plot_stems[collisions] <- paste0(plot_stems[collisions], "__", collision_hashes)

    for (index in seq_len(nrow(res))) {
        plot_path <- file.path(plots_dir, paste0(plot_stems[[index]], ".png"))
        device_open <- FALSE
        plot_error <- tryCatch(
            {
                png(plot_path, width = 1200, height = 800, res = 130)
                device_open <- TRUE
                hist(
                    reference_values,
                    breaks = "FD",
                    col = "#8CB9D9",
                    border = "white",
                    main = paste("PGS reference distribution:", res$IID[[index]]),
                    xlab = "Raw PGS (PLINK --score sum)"
                )
                abline(v = res$PGS[[index]], col = "#C43C39", lwd = 3)
                legend(
                    "topright",
                    legend = sprintf(
                        "Sample PGS %.5f | percentile %.2f",
                        res$PGS[[index]],
                        res$PERCENTILE[[index]]
                    ),
                    col = "#C43C39",
                    lwd = 3,
                    bty = "n"
                )
                if (res$QC_STATUS[[index]] == "WARN") {
                    mtext(
                        paste("WARNING:", res$WARNING_CODE[[index]]),
                        side = 3,
                        col = "#B22222",
                        line = 0.2,
                        cex = 0.75
                    )
                }
                mtext(
                    "Research use only. No clinical risk interpretation is provided.",
                    side = 1,
                    line = 3,
                    cex = 0.75
                )
                dev.off()
                device_open <- FALSE
                NULL
            },
            error = function(error) error
        )
        if (inherits(plot_error, "error")) {
            if (device_open) {
                try(dev.off(), silent = TRUE)
            }
            unlink(plot_path, force = TRUE)
            prior_code <- res$WARNING_CODE[[index]]
            res$QC_STATUS[[index]] <- "WARN"
            res$WARNING_CODE[[index]] <- if (nzchar(prior_code)) {
                paste0(prior_code, ";PLOT_FAILED")
            } else {
                "PLOT_FAILED"
            }
            res$WARNING_MESSAGE[[index]] <- paste(
                res$WARNING_MESSAGE[[index]],
                warning_text[["PLOT_FAILED"]]
            )
        }
    }
}
write_atomic(res, results_path)

# Assemble the English run report, including checksums, QC summaries, commands, and PLINK logs.
palindromic_count <- model_metadata |>
    filter(.data$field == "palindromic_nonzero_rows") |>
    pull("value")
warning_codes <- if (all(res$QC_STATUS == "PASS")) {
    "None"
} else {
    unique(unlist(strsplit(res$WARNING_CODE[nzchar(res$WARNING_CODE)], ";", fixed = TRUE)))
}
report <- c(
    paste0("LN-PriGI-PRS run report | version ", VERSION),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "RESEARCH-USE NOTICE",
    paste(
        "This software reports PGS percentiles relative to a specified reference population.",
        "It does not provide clinical risk interpretation, diagnosis, or treatment advice."
    ),
    "",
    "RUN CONFIGURATION",
    paste0("Input prefix: ", input_prefix),
    paste0("Input format: ", input_mode),
    "Genome build: hg19",
    paste0("Output directory: ", out_dir),
    paste0("Output prefix: ", out_prefix),
    paste0("Plot mode: ", plot_mode),
    paste0("Overwrite: ", isTRUE(opt$overwrite)),
    paste0("PLINK path: ", plink),
    paste0("PLINK version: ", paste(version_lines, collapse = " ")),
    paste0("Reference: ", reference$source),
    paste0("Reference n: ", reference_n),
    sprintf("Reference adjusted Fisher-Pearson skewness G1: %.6f", reference$skew),
    paste0("MAF filter: ", if (is.null(maf)) "not used" else maf),
    paste0("GENO filter: ", if (is.null(geno)) "not used" else geno),
    "",
    "FROZEN ASSETS",
    paste0("Model version: ", VERSION),
    paste0(manifest$path, " SHA-256: ", manifest$sha256),
    "",
    "MODEL AND ALIGNMENT",
    paste0("Eligible nonzero-weight model variants: ", model_n),
    paste0("Input duplicate coordinates excluded: ", nrow(dup_coord)),
    paste0("Input duplicate variant IDs excluded: ", nrow(dup_id)),
    paste0("Direct allele matches: ", direct_match_n),
    paste0("Complement allele matches: ", complement_match_n),
    paste0("Allele-incompatible or unmatched input variants: ", nrow(bim) - nrow(matched)),
    paste0("Permanently excluded palindromic nonzero model variants: ", palindromic_count),
    paste0("Variants retained for scoring: ", retained_n),
    sprintf("Retained absolute-weight coverage: %.6f", abs_coverage),
    "",
    "WARNINGS",
    warning_codes,
    "",
    "COMMANDS",
    commands,
    logs,
    "",
    "OUTPUTS",
    paste0("Results: ", results_path),
    paste0("Plots: ", if (plot_mode == "all") plots_dir else "disabled"),
    "Completion status: SUCCESS"
)
report_tmp <- paste0(report_path, ".tmp")
writeLines(report, report_tmp, useBytes = TRUE)
if (!file.rename(report_tmp, report_path)) {
    fail("Cannot finalize report.")
}
cat("Completed. Results: ", results_path, "\n", sep = "")
