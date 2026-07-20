<#
.SYNOPSIS
Packages the fixed PLINK 1 acceptance genotypes and refreshes SHA-256 manifests.

.DESCRIPTION
Maintainer-only release tool. Routine users do not need to run it. Supply the
prefix of the fixed BED/BIM/FAM test cohort after rebuilding frozen assets.

.EXAMPLE
pwsh tools/package_fixed_test_genotypes.ps1 -TestPrefix D:\data\dummy_genome
#>
param(
    [Parameter(Mandatory = $true)][string]$TestPrefix,
    [string]$DataDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) "data")
)

$ErrorActionPreference = "Stop"
$testDirectory = Join-Path $DataDirectory "test_data"
$archive = Join-Path $testDirectory "dummy_genome_plink1.zip"
$files = @(".bed", ".bim", ".fam") | ForEach-Object { "$TestPrefix$_" }
foreach ($file in $files) {
    if (-not (Test-Path -LiteralPath $file -PathType Leaf)) { throw "Missing PLINK test file: $file" }
}
# Package the unchanged acceptance genotypes into the repository release asset.
Compress-Archive -LiteralPath $files -DestinationPath $archive -CompressionLevel Optimal -Force
$archiveRows = @("sha256`tpath")
foreach ($file in $files) {
    $hash = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant()
    $archiveRows += "$hash`t$([IO.Path]::GetFileName($file))"
}
$archiveManifest = Join-Path $testDirectory "archive_contents.sha256.tsv"
[IO.File]::WriteAllText($archiveManifest, ([string]::Join([char]10, $archiveRows) + [char]10), [Text.UTF8Encoding]::new($false))

# Refresh the top-level manifest only after the archive-content hashes exist.
$relativeAssets = @(
    "model/gwas_prior_union.model.rds",
    "model/model_metadata.tsv",
    "reference/reference_pgs.tsv",
    "reference/reference_metadata.tsv",
    "test_data/dummy_genome_plink1.zip",
    "test_data/archive_contents.sha256.tsv",
    "test_data/expected_pgs.tsv"
)
$rows = @("sha256`tpath")
foreach ($relative in $relativeAssets) {
    $path = Join-Path $DataDirectory ($relative.Replace("/", [IO.Path]::DirectorySeparatorChar))
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Missing frozen asset: $relative" }
    $hash = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
    $rows += "$hash`t$relative"
}
$manifest = Join-Path $DataDirectory "manifest.sha256.tsv"
[IO.File]::WriteAllText($manifest, ([string]::Join([char]10, $rows) + [char]10), [Text.UTF8Encoding]::new($false))
Write-Output "Packaged and checksummed: $archive"
