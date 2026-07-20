#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(digest)})

VERSION <- "0.1.0"
args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
ROOT <- normalizePath(file.path(dirname(script_arg), ".."), winslash = "/", mustWork = TRUE)

fail <- function(...) stop(paste0(...), call. = FALSE)
"%||%" <- function(x, y) if (is.null(x)) y else x
parse_cli <- function(x) {
  flags <- c("--overwrite", "--help", "--version")
  out <- list(); i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) fail("Unexpected argument: ", key)
    if (key %in% flags) { out[[substring(key, 3L)]] <- TRUE; i <- i + 1L; next }
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) fail("Missing value for ", key)
    out[[substring(key, 3L)]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
usage <- function() cat(paste0(
"LN-PriGI-PRS ", VERSION, "\n\n",
"Usage:\n  Rscript R/main.R --input PREFIX --genome-build hg19 --output-dir DIR [options]\n\n",
"Required:\n  --input PATH       PLINK 1 BED/BIM/FAM or PED/MAP prefix\n",
"  --genome-build hg19      Explicit genome-build declaration; only hg19 is accepted\n",
"  --output-dir DIR         Output directory\n\n",
"Options:\n  --output-prefix NAME      Output prefix [ln_prs]\n",
"  --reference-pgs FILE     Replacement reference TSV: sample_id, pgs_score\n",
"  --plink PATH             PLINK 1.9 executable [search PATH]\n",
"  --maf FLOAT              Optional PLINK MAF filter; no default\n",
"  --geno FLOAT             Optional PLINK missingness filter; no default\n",
"  --plot-mode all|none     Per-sample plots [all]\n",
"  --overwrite              Replace this exact output prefix\n",
"  --help                    Show help\n  --version                 Show version\n"))

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) { usage(); quit(status = 0L) }
if (isTRUE(opt$version)) { cat(VERSION, "\n"); quit(status = 0L) }
if ("input-prefix" %in% names(opt) && !"input" %in% names(opt)) opt$input <- opt[["input-prefix"]]
required <- c("input", "genome-build", "output-dir")
missing_opt <- required[!required %in% names(opt)]
if (length(missing_opt)) fail("Missing required option(s): --", paste(missing_opt, collapse = ", --"))
if (!identical(opt[["genome-build"]], "hg19")) fail("Only hg19 input is accepted. Liftover is not provided.")
plot_mode <- opt[["plot-mode"]] %||% "all"
if (!plot_mode %in% c("all", "none")) fail("--plot-mode must be 'all' or 'none'.")
out_prefix <- opt[["output-prefix"]] %||% "ln_prs"
if (!grepl("^[A-Za-z0-9._-]+$", out_prefix)) fail("--output-prefix may contain only letters, digits, dot, underscore, and hyphen.")
validate_num <- function(name, lower, upper) {
  if (!name %in% names(opt)) return(NULL)
  z <- suppressWarnings(as.numeric(opt[[name]]))
  if (!is.finite(z) || z < lower || z > upper) fail("--", name, " must be between ", lower, " and ", upper, ".")
  z
}
maf <- validate_num("maf", 0, 0.5); geno <- validate_num("geno", 0, 1)

input_prefix <- normalizePath(opt[["input"]], winslash = "/", mustWork = FALSE)
has <- function(ext) file.exists(paste0(input_prefix, ext))
binary <- all(vapply(c(".bed", ".bim", ".fam"), has, logical(1)))
text <- all(vapply(c(".ped", ".map"), has, logical(1)))
any_binary <- any(vapply(c(".bed", ".bim", ".fam"), has, logical(1)))
any_text <- any(vapply(c(".ped", ".map"), has, logical(1)))
if (binary && text) fail("Both BED/BIM/FAM and PED/MAP were found. Provide an unambiguous prefix.")
if (!binary && !text) fail(if (any_binary || any_text) "The PLINK input file set is incomplete." else "No supported PLINK 1 input files were found.")
input_mode <- if (binary) "BED/BIM/FAM" else "PED/MAP"

plink <- opt$plink %||% Sys.which("plink")
if (!nzchar(plink) || !file.exists(plink)) fail("PLINK was not found. Install PLINK 1.9 or provide --plink PATH.")
plink <- normalizePath(plink, winslash = "/", mustWork = TRUE)
version_lines <- suppressWarnings(system2(plink, "--version", stdout = TRUE, stderr = TRUE))
if (!any(grepl("PLINK v1\\.9", version_lines))) fail("The selected executable is not PLINK 1.9: ", paste(version_lines, collapse = " "))

out_dir <- normalizePath(opt[["output-dir"]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) fail("Cannot create output directory: ", out_dir)
results_path <- file.path(out_dir, paste0(out_prefix, ".results.tsv"))
report_path <- file.path(out_dir, paste0(out_prefix, ".run_report.txt"))
plots_dir <- file.path(out_dir, paste0(out_prefix, ".plots"))
targets <- c(results_path, report_path, plots_dir)
existing <- targets[file.exists(targets) | dir.exists(targets)]
if (length(existing) && !isTRUE(opt$overwrite)) fail("Output already exists. Use --overwrite to replace this exact prefix: ", existing[[1]])
if (isTRUE(opt$overwrite)) {
  for (z in existing) if (dir.exists(z)) unlink(z, recursive = TRUE, force = TRUE) else unlink(z, force = TRUE)
}
work <- tempfile(paste0("ln_prs_", out_prefix, "_"), tmpdir = out_dir)
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
logs <- character(); commands <- character()
run_plink <- function(label, x) {
  log <- file.path(work, paste0(label, ".log")); commands <<- c(commands, paste(shQuote(plink), paste(shQuote(x), collapse = " ")))
  status <- system2(plink, x, stdout = log, stderr = log)
  txt <- readLines(log, warn = FALSE); logs <<- c(logs, paste0("\n===== PLINK: ", label, " =====\n"), txt)
  if (!identical(status, 0L)) {
    writeLines(c("LN-PriGI-PRS failed.", paste0("Stage: ", label), "", commands, logs), report_path, useBytes = TRUE)
    fail("PLINK failed during ", label, ". Details were written to ", report_path)
  }
  txt
}

manifest_path <- file.path(ROOT, "data/manifest.sha256.tsv")
manifest <- fread(manifest_path)
required_assets <- c("model/gwas_prior_union.model.rds", "model/model_metadata.tsv")
if (!"reference-pgs" %in% names(opt)) required_assets <- c(required_assets, "reference/reference_pgs.tsv", "reference/reference_metadata.tsv")
manifest <- manifest[path %chin% required_assets]
for (i in seq_len(nrow(manifest))) {
  asset <- file.path(ROOT, "data", manifest$path[[i]])
  if (!file.exists(asset) || digest(asset, algo = "sha256", file = TRUE) != manifest$sha256[[i]]) fail("Asset checksum verification failed: ", manifest$path[[i]])
}
model_metadata <- fread(file.path(ROOT, "data/model/model_metadata.tsv"))
model <- readRDS(file.path(ROOT, "data/model/gwas_prior_union.model.rds"))
setDT(model); model[, `:=`(EFFECT_ALLELE = toupper(EFFECT_ALLELE), OTHER_ALLELE = toupper(OTHER_ALLELE))]
model_n <- nrow(model); model_abs <- sum(abs(model$WEIGHT))

base_prefix <- input_prefix
if (text) {
  converted <- file.path(work, "ped_converted")
  run_plink("convert_ped", c("--file", input_prefix, "--make-bed", "--out", converted))
  base_prefix <- converted
}
bim <- fread(paste0(base_prefix, ".bim"), col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2"))
bim[, `:=`(A1 = toupper(A1), A2 = toupper(A2), POS = as.integer(POS))]
dup_coord <- bim[, .N, by = .(CHR, POS)][N > 1L, .(CHR, POS)]
if (nrow(dup_coord)) bim <- bim[!dup_coord, on = .(CHR, POS)]
dup_id <- bim[, .N, by = SNP_ID][N > 1L, .(SNP_ID)]
if (nrow(dup_id)) bim <- bim[!dup_id, on = "SNP_ID"]
matched <- merge(bim[, .(SNP_ID, CHR, POS, USER_A1 = A1, USER_A2 = A2)], model, by = c("CHR", "POS"), sort = FALSE)
comp <- function(x) chartr("ACGT", "TGCA", x)
matched[, DIRECT_MATCH := (EFFECT_ALLELE == USER_A1 & OTHER_ALLELE == USER_A2) | (EFFECT_ALLELE == USER_A2 & OTHER_ALLELE == USER_A1)]
matched[, COMPLEMENT_MATCH := (comp(EFFECT_ALLELE) == USER_A1 & comp(OTHER_ALLELE) == USER_A2) | (comp(EFFECT_ALLELE) == USER_A2 & comp(OTHER_ALLELE) == USER_A1)]
matched <- matched[DIRECT_MATCH | COMPLEMENT_MATCH]
matched[, USER_EFFECT := fifelse(DIRECT_MATCH, EFFECT_ALLELE, comp(EFFECT_ALLELE))]
matched[, USER_OTHER := fifelse(USER_EFFECT == USER_A1, USER_A2, USER_A1)]
if (!nrow(matched)) fail("No eligible model variants matched the input genotype data.")
direct_match_n <- matched[DIRECT_MATCH == TRUE, .N]
complement_match_n <- matched[COMPLEMENT_MATCH == TRUE & DIRECT_MATCH == FALSE, .N]
score0 <- file.path(work, "score_initial.tsv"); fwrite(matched[, .(SNP_ID, USER_EFFECT, WEIGHT)], score0, sep = "\t", quote = FALSE)
qc <- file.path(work, "qc")
qc_args <- c("--bfile", base_prefix, "--extract", score0)
if (!is.null(maf)) qc_args <- c(qc_args, "--maf", format(maf, scientific = FALSE))
if (!is.null(geno)) qc_args <- c(qc_args, "--geno", format(geno, scientific = FALSE))
run_plink("extract_and_filter", c(qc_args, "--make-bed", "--out", qc))
qc_bim <- fread(paste0(qc, ".bim"), col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2"))
retained <- matched[SNP_ID %chin% qc_bim$SNP_ID]
if (!nrow(retained)) fail("No model variants remained after optional filters.")
score <- file.path(work, "score.tsv"); a2 <- file.path(work, "a2.tsv")
fwrite(retained[, .(SNP_ID, USER_EFFECT, WEIGHT)], score, sep = "\t", quote = FALSE)
fwrite(retained[, .(SNP_ID, USER_OTHER)], a2, sep = "\t", quote = FALSE, col.names = FALSE)
missing_prefix <- file.path(work, "missing")
run_plink("missingness", c("--bfile", qc, "--missing", "--out", missing_prefix))
imiss <- fread(paste0(missing_prefix, ".imiss"))
if (anyDuplicated(imiss[, .(FID, IID)])) fail("FID/IID combinations must be unique.")
retained_n <- nrow(retained)
imiss[, NONZERO_SNPS_USED := pmax(0, retained_n - N_MISS)]
imiss[, NONZERO_VARIANT_COVERAGE := NONZERO_SNPS_USED / model_n]
abs_coverage <- sum(abs(retained$WEIGHT)) / model_abs
score_ready <- file.path(work, "score_ready")
a2_log <- run_plink("encode_missing_as_zero", c("--bfile", qc, "--a2-allele", a2, "2", "1", "--fill-missing-a2", "--make-bed", "--out", score_ready))
if (any(grepl("Impossible A2 allele assignment", a2_log, fixed = TRUE))) fail("Internal allele alignment verification failed during A2 assignment.")
prs_prefix <- file.path(work, "prs")
run_plink("score", c("--bfile", score_ready, "--score", score, "1", "2", "3", "header", "sum", "--out", prs_prefix))
prs <- fread(paste0(prs_prefix, ".profile"))[, .(FID, IID, PGS = SCORESUM)]
res <- merge(prs, imiss[, .(FID, IID, NONZERO_SNPS_USED, NONZERO_VARIANT_COVERAGE)], by = c("FID", "IID"), sort = FALSE)
res[FID == "0", FID := NA_character_]

read_reference <- function(path = NULL) {
  src <- if (is.null(path)) file.path(ROOT, "data/reference/reference_pgs.tsv") else path
  x <- fread(src, sep = "\t")
  if (!identical(names(x), c("sample_id", "pgs_score"))) fail("Reference PGS must be a tab-delimited table with exactly: sample_id, pgs_score.")
  if (nrow(x) < 100L) fail("The reference PGS has too few samples to represent population distribution characteristics (n < 100).")
  if (anyNA(x$sample_id) || any(!nzchar(trimws(x$sample_id))) || anyDuplicated(x$sample_id)) fail("Reference sample_id values must be nonblank and unique.")
  x[, pgs_score := suppressWarnings(as.numeric(pgs_score))]
  if (any(!is.finite(x$pgs_score)) || sd(x$pgs_score) == 0) fail("Reference pgs_score values must be finite and have nonzero variance.")
  n <- nrow(x); z <- (x$pgs_score - mean(x$pgs_score)) / sd(x$pgs_score)
  skew <- n / ((n - 1) * (n - 2)) * sum(z^3)
  list(data = x, source = if (is.null(path)) "built-in:LN_WGS_controls" else normalizePath(path, winslash = "/", mustWork = TRUE), skew = skew)
}
ref <- read_reference(opt[["reference-pgs"]]); rv <- ref$data$pgs_score; rn <- length(rv)
res[, PERCENTILE := round(vapply(PGS, function(v) 100 * (sum(rv < v) + 0.5 * sum(rv == v)) / rn, numeric(1)), 2)]
res[, `:=`(NONZERO_MODEL_SNPS = model_n, ABS_WEIGHT_COVERAGE = abs_coverage)]
res[, warnings__ := ""]
add_warning <- function(code, idx = rep(TRUE, nrow(res))) res[idx, warnings__ := fifelse(nzchar(warnings__), paste0(warnings__, ";", code), code)]
if (rn < 200L) add_warning("REFERENCE_SAMPLE_SIZE_LOW")
if (abs(ref$skew) > 1) add_warning("REFERENCE_DISTRIBUTION_SKEWED")
add_warning("NONZERO_VARIANT_COVERAGE_BELOW_90_PERCENT", res$NONZERO_VARIANT_COVERAGE < 0.90)
warning_text <- c(
  REFERENCE_SAMPLE_SIZE_LOW = "The reference sample size is below 200 and may not fully represent the population distribution.",
  REFERENCE_DISTRIBUTION_SKEWED = "The reference PGS distribution has absolute adjusted Fisher-Pearson skewness above 1.",
  NONZERO_VARIANT_COVERAGE_BELOW_90_PERCENT = "Nonzero-weight model variant coverage is below 90%.",
  PLOT_FAILED = "The per-sample plot could not be generated."
)
res[, `:=`(GENOME_BUILD = "hg19", REFERENCE_SOURCE = ref$source,
           PERCENTILE_RELIABLE = NONZERO_VARIANT_COVERAGE >= 0.90,
           QC_STATUS = fifelse(nzchar(warnings__), "WARN", "PASS"),
           WARNING_CODE = warnings__)]
res[, WARNING_MESSAGE := vapply(strsplit(WARNING_CODE, ";", fixed = TRUE), function(z) paste(unname(warning_text[z[nzchar(z)]]), collapse = " "), character(1))]
res[, warnings__ := NULL]
setcolorder(res, c("FID", "IID", "PGS", "PERCENTILE", "NONZERO_MODEL_SNPS", "NONZERO_SNPS_USED", "NONZERO_VARIANT_COVERAGE", "ABS_WEIGHT_COVERAGE", "QC_STATUS", "PERCENTILE_RELIABLE", "WARNING_CODE", "WARNING_MESSAGE", "REFERENCE_SOURCE", "GENOME_BUILD"))

write_atomic <- function(x, path) { tmp <- paste0(path, ".tmp"); fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA"); if (!file.rename(tmp, path)) fail("Cannot finalize output: ", path) }
if (plot_mode == "all") {
  dir.create(plots_dir)
  safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
  plot_stems <- paste0(safe(res$FID), "__", safe(res$IID))
  collisions <- duplicated(plot_stems) | duplicated(plot_stems, fromLast = TRUE)
  collision_hashes <- substr(vapply(paste(res$FID[collisions], res$IID[collisions], sep = "\t"), digest, character(1), algo = "sha256", serialize = FALSE), 1L, 8L)
  plot_stems[collisions] <- paste0(plot_stems[collisions], "__", collision_hashes)
  for (i in seq_len(nrow(res))) {
    f <- file.path(plots_dir, paste0(plot_stems[[i]], ".png"))
    device_open <- FALSE
    plot_error <- tryCatch({
      png(f, width = 1200, height = 800, res = 130); device_open <- TRUE
      hist(rv, breaks = "FD", col = "#8CB9D9", border = "white", main = paste("PGS reference distribution:", res$IID[[i]]), xlab = "Raw PGS (PLINK --score sum)")
      abline(v = res$PGS[[i]], col = "#C43C39", lwd = 3)
      legend("topright", legend = sprintf("Sample PGS %.5f | percentile %.2f", res$PGS[[i]], res$PERCENTILE[[i]]), col = "#C43C39", lwd = 3, bty = "n")
      if (res$QC_STATUS[[i]] == "WARN") mtext(paste("WARNING:", res$WARNING_CODE[[i]]), side = 3, col = "#B22222", line = 0.2, cex = 0.75)
      mtext("Research use only. No clinical risk interpretation is provided.", side = 1, line = 3, cex = 0.75)
      dev.off(); device_open <- FALSE
      NULL
    }, error = function(e) e)
    if (inherits(plot_error, "error")) {
      if (device_open) try(dev.off(), silent = TRUE)
      unlink(f, force = TRUE)
      res[i, `:=`(QC_STATUS = "WARN", WARNING_CODE = fifelse(nzchar(WARNING_CODE), paste0(WARNING_CODE, ";PLOT_FAILED"), "PLOT_FAILED"), WARNING_MESSAGE = paste(WARNING_MESSAGE, warning_text[["PLOT_FAILED"]]))]
    }
  }
}
write_atomic(res, results_path)
report <- c(
  paste0("LN-PriGI-PRS run report | version ", VERSION),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "", "RESEARCH-USE NOTICE", "This software reports PGS percentiles relative to a specified reference population. It does not provide clinical risk interpretation, diagnosis, or treatment advice.",
  "", "RUN CONFIGURATION", paste0("Input prefix: ", input_prefix), paste0("Input format: ", input_mode), "Genome build: hg19", paste0("Output directory: ", out_dir), paste0("Output prefix: ", out_prefix), paste0("Plot mode: ", plot_mode), paste0("Overwrite: ", isTRUE(opt$overwrite)), paste0("PLINK path: ", plink), paste0("PLINK version: ", paste(version_lines, collapse = " ")), paste0("Reference: ", ref$source), paste0("Reference n: ", rn), sprintf("Reference adjusted Fisher-Pearson skewness G1: %.6f", ref$skew), paste0("MAF filter: ", if (is.null(maf)) "not used" else maf), paste0("GENO filter: ", if (is.null(geno)) "not used" else geno),
  "", "FROZEN ASSETS", paste0("Model version: ", VERSION), paste0(manifest$path, " SHA-256: ", manifest$sha256),
  "", "MODEL AND ALIGNMENT", paste0("Eligible nonzero-weight model variants: ", model_n), paste0("Input duplicate coordinates excluded: ", nrow(dup_coord)), paste0("Input duplicate variant IDs excluded: ", nrow(dup_id)), paste0("Direct allele matches: ", direct_match_n), paste0("Complement allele matches: ", complement_match_n), paste0("Allele-incompatible or unmatched input variants: ", nrow(bim) - nrow(matched)), paste0("Permanently excluded palindromic nonzero model variants: ", model_metadata[field == "palindromic_nonzero_rows", value]), paste0("Variants retained for scoring: ", retained_n), sprintf("Retained absolute-weight coverage: %.6f", abs_coverage),
  "", "WARNINGS", if (all(res$QC_STATUS == "PASS")) "None" else unique(unlist(strsplit(res$WARNING_CODE[nzchar(res$WARNING_CODE)], ";", fixed = TRUE))),
  "", "COMMANDS", commands, logs,
  "", "OUTPUTS", paste0("Results: ", results_path), paste0("Plots: ", if (plot_mode == "all") plots_dir else "disabled"), "Completion status: SUCCESS"
)
report_tmp <- paste0(report_path, ".tmp"); writeLines(report, report_tmp, useBytes = TRUE); if (!file.rename(report_tmp, report_path)) fail("Cannot finalize report.")
cat("Completed. Results: ", results_path, "\n", sep = "")
