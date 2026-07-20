#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: build_assets.R <weight> <harmonized_bim> <reference_profile> <test_prefix> <output_data_dir> <plink>", call. = FALSE)
}
weight_path <- normalizePath(args[[1]], mustWork = TRUE)
bim_path <- normalizePath(args[[2]], mustWork = TRUE)
reference_path <- normalizePath(args[[3]], mustWork = TRUE)
test_prefix <- normalizePath(args[[4]], mustWork = FALSE)
out_dir <- normalizePath(args[[5]], mustWork = FALSE)
plink <- normalizePath(args[[6]], mustWork = TRUE)

dir.create(file.path(out_dir, "model"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "reference"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "test_data"), recursive = TRUE, showWarnings = FALSE)

complement <- function(x) chartr("ATGC", "TACG", x)
is_dna <- function(x) !is.na(x) & nchar(x) == 1L & x %chin% c("A", "C", "G", "T")
is_palindromic <- function(a1, a2) paste0(a1, a2) %chin% c("AT", "TA", "CG", "GC")

message("Reading model weights...")
weights <- fread(weight_path, select = c("posID", "A1", "A1Effect"))
weights[, c("CHR", "POS") := tstrsplit(sub("\\[hg19\\]$", "", posID), ":", fixed = TRUE)]
weights[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS), A1 = toupper(A1), A1Effect = as.numeric(A1Effect))]
weights[, COORD := paste(CHR, POS, sep = ":")]
weight_dup_coords <- weights[duplicated(COORD) | duplicated(COORD, fromLast = TRUE), unique(COORD)]
weights_unique <- weights[!COORD %chin% weight_dup_coords]

message("Reading harmonized BIM allele source...")
bim <- fread(bim_path, col.names = c("CHR", "SNP_ID", "CM", "POS", "B1", "B2"))[, .(CHR, SNP_ID, POS, B1, B2)]
bim[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS), B1 = toupper(B1), B2 = toupper(B2))]
bim[, COORD := paste(CHR, POS, sep = ":")]
bim_dup_coords <- bim[duplicated(COORD) | duplicated(COORD, fromLast = TRUE), unique(COORD)]
bim_unique <- bim[!COORD %chin% bim_dup_coords]

message("Recovering the second model allele...")
joined <- merge(weights_unique, bim_unique, by = c("CHR", "POS", "COORD"), all.x = TRUE, sort = FALSE)
joined[, DIRECT := is_dna(A1) & is_dna(B1) & is_dna(B2) & (A1 == B1 | A1 == B2)]
joined[, COMPLEMENT := !DIRECT & is_dna(A1) & is_dna(B1) & is_dna(B2) & (complement(A1) == B1 | complement(A1) == B2)]
joined[, A2 := fifelse(DIRECT & A1 == B1, B2,
                fifelse(DIRECT & A1 == B2, B1,
                fifelse(COMPLEMENT & complement(A1) == B1, complement(B2),
                fifelse(COMPLEMENT & complement(A1) == B2, complement(B1), NA_character_))))]
joined[, VALID_PAIR := is_dna(A1) & is_dna(A2) & A1 != A2]
joined[, PALINDROMIC := VALID_PAIR & is_palindromic(A1, A2)]
joined[, NONZERO := is.finite(A1Effect) & A1Effect != 0]

model <- joined[VALID_PAIR & NONZERO & !PALINDROMIC,
  .(MODEL_ID = posID, CHR, POS, EFFECT_ALLELE = A1, OTHER_ALLELE = A2, WEIGHT = A1Effect)]
setorder(model, CHR, POS, MODEL_ID)
if (nrow(model) == 0L) stop("No eligible nonzero model variants were built.", call. = FALSE)

model_path <- file.path(out_dir, "model", "gwas_prior_union.model.rds")
saveRDS(model, model_path, compress = "xz", version = 3)
metadata <- data.table(
  field = c("model_name", "genome_build", "source_weight_file", "allele_source_file",
          "weight_rows", "zero_weight_rows", "nonzero_weight_rows", "duplicate_weight_coordinates",
          "duplicate_allele_source_coordinates", "unmatched_or_invalid_rows", "palindromic_nonzero_rows",
          "eligible_nonzero_rows", "model_rds_sha256"),
  value = c("gwas_prior_union", "hg19", basename(weight_path), basename(bim_path),
            nrow(weights), weights[A1Effect == 0, .N], weights[is.finite(A1Effect) & A1Effect != 0, .N],
            length(weight_dup_coords), length(bim_dup_coords), joined[VALID_PAIR == FALSE, .N],
            joined[VALID_PAIR & NONZERO & PALINDROMIC, .N], nrow(model),
            digest(model_path, algo = "sha256", file = TRUE))
)
fwrite(metadata, file.path(out_dir, "model", "model_metadata.tsv"), sep = "\t", quote = FALSE)

message("Building the built-in control reference...")
reference_profile <- fread(reference_path)
required_reference <- c("FID", "IID", "SCORESUM")
if (!all(required_reference %chin% names(reference_profile))) stop("Reference profile is missing FID/IID/SCORESUM.", call. = FALSE)
reference <- reference_profile[grepl("^CTRL", FID), .(sample_id = as.character(IID), pgs_score = as.numeric(SCORESUM))]
if (nrow(reference) == 0L) stop("No CTRL reference samples found.", call. = FALSE)
reference_out <- file.path(out_dir, "reference", "reference_pgs.tsv")
fwrite(reference, reference_out, sep = "\t", quote = FALSE)
reference_meta <- data.table(
  field = c("reference_name", "source_profile", "selection", "sample_count", "reference_sha256"),
  value = c("LN_WGS_controls", basename(reference_path), "FID starts with CTRL", nrow(reference),
            digest(reference_out, algo = "sha256", file = TRUE))
)
fwrite(reference_meta, file.path(out_dir, "reference", "reference_metadata.tsv"), sep = "\t", quote = FALSE)

message("Preparing score files for the fixed test genotypes...")
test_bim <- fread(paste0(test_prefix, ".bim"), col.names = c("CHR", "SNP_ID", "CM", "POS", "A1", "A2"))
test_bim[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
test_matches <- merge(test_bim[, .(SNP_ID, CHR, POS, TEST_A1 = A1, TEST_A2 = A2)], model,
                      by = c("CHR", "POS"), sort = FALSE)
test_matches[, USER_EFFECT := fifelse(EFFECT_ALLELE == TEST_A1 | EFFECT_ALLELE == TEST_A2, EFFECT_ALLELE, complement(EFFECT_ALLELE))]
test_matches <- test_matches[USER_EFFECT == TEST_A1 | USER_EFFECT == TEST_A2]
test_matches[, USER_OTHER := fifelse(USER_EFFECT == TEST_A1, TEST_A2, TEST_A1)]
score_file <- file.path(out_dir, "test_data", "test_score.tsv")
a2_file <- file.path(out_dir, "test_data", "test_a2.tsv")
fwrite(test_matches[, .(SNP_ID, USER_EFFECT, WEIGHT)], score_file, sep = "\t", quote = FALSE)
fwrite(test_matches[, .(SNP_ID, USER_OTHER)], a2_file, sep = "\t", quote = FALSE, col.names = FALSE)

work <- tempfile("ln_prs_asset_build_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
plink_log <- file.path(work, "plink_asset_build.log")
run_plink <- function(args) {
  status <- system2(plink, args, stdout = plink_log, stderr = plink_log)
  if (!identical(status, 0L)) stop(paste(readLines(plink_log, warn = FALSE), collapse = "\n"), call. = FALSE)
}
extracted <- file.path(work, "extracted")
run_plink(c("--bfile", test_prefix, "--extract", score_file, "--make-bed", "--out", extracted))
score_ready <- file.path(work, "score_ready")
run_plink(c("--bfile", extracted, "--a2-allele", a2_file, "2", "1", "--fill-missing-a2",
            "--make-bed", "--out", score_ready))
expected_prefix <- file.path(work, "expected")
run_plink(c("--bfile", score_ready, "--score", score_file, "1", "2", "3", "header", "sum",
            "--out", expected_prefix))
expected <- fread(paste0(expected_prefix, ".profile"))[, .(FID, IID, EXPECTED_PGS = SCORESUM)]
fwrite(expected, file.path(out_dir, "test_data", "expected_pgs.tsv"), sep = "\t", quote = FALSE)

message("Writing SHA-256 manifests...")
assets <- c(
  file.path("model", "gwas_prior_union.model.rds"),
  file.path("model", "model_metadata.tsv"),
  file.path("reference", "reference_pgs.tsv"),
  file.path("reference", "reference_metadata.tsv"),
  file.path("test_data", "expected_pgs.tsv")
)
manifest <- rbindlist(lapply(assets, function(rel) {
  path <- file.path(out_dir, rel)
  data.table(sha256 = digest(path, algo = "sha256", file = TRUE), path = gsub("\\\\", "/", rel))
}))
fwrite(manifest, file.path(out_dir, "manifest.sha256.tsv"), sep = "\t", quote = FALSE)
message("Asset build complete: ", out_dir)
