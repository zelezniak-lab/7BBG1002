#!/bin/bash -l
#SBATCH --job-name=install-rpkgs
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

# ---- Modules (same stack used in teaching sessions) ----
module load rstudio/v2023.03.0_386-gcc-13.2.0-r-4.3.0-python-3.11.6
module load freetype/2.11.1-gcc-13.2.0
module load libpng/1.6.39-gcc-13.2.0-curl-8.4.0
module load libtiff/4.5.1-gcc-13.2.0-curl-8.4.0
module load libwebp/1.2.4-gcc-13.2.0
module load libjpeg-turbo/3.0.0-gcc-13.2.0-curl-8.4.0
module load pkgconf/1.9.5-gcc-13.2.0
module load cmake/3.27.7-gcc-13.2.0
module load libxml2/2.10.3-gcc-13.2.0
module load openblas/0.3.24-gcc-13.2.0
module load intel-mkl/2020.4.304-gcc-13.2.0-openmp
module load netlib-scalapack/2.2.0-gcc-13.2.0-openmpi-4.1.6-python-3.11.6
module load openssl/3.1.3-gcc-13.2.0
module load fontconfig/2.14.2-gcc-13.2.0-python-3.11.6
# (If builds complain about curl: module load curl)

export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
export PKG_CONFIG_ALLOW_SYSTEM_LIBS=1

# ---- Runtime linker paths via pkg-config (no EasyBuild) ----
for pc in freetype2 libpng libtiff-4 libjpeg libwebp libwebpmux; do
  d=$(pkg-config --variable=libdir "$pc" 2>/dev/null || true)
  [[ -n "$d" ]] && export LD_LIBRARY_PATH="$d:$LD_LIBRARY_PATH"
done

# ---- Shared site library ----
export R_VERSION_SHORT=4.3
export R_LIBS_SITE="/scratch/grp/msc_appbio/shared/Rlib/${R_VERSION_SHORT}"
mkdir -p "$R_LIBS_SITE"

echo "Installing R + Bioconductor packages into: $R_LIBS_SITE"
module list

# ---- Unlock (chmod only) ----
echo "🔓 Unlocking $R_LIBS_SITE..."
chmod -R u+rwX,g+rwX "$R_LIBS_SITE" 2>/dev/null || true
umask 002

# ---- Run R non-interactively ----
Rscript --vanilla - <<'RSCRIPT'
site_lib <- Sys.getenv("R_LIBS_SITE")
repos_cran <- "https://cran.ma.imperial.ac.uk"

# Put site lib first so everything installs/loads there
.libPaths(c(site_lib, .libPaths()))
options(repos = c(CRAN = repos_cran))

message("Site library: ", site_lib)
dir.create(site_lib, recursive = TRUE, showWarnings = FALSE)

# ----- Ensure BiocManager in site lib -----
if (!requireNamespace("BiocManager", quietly = TRUE, lib.loc = site_lib)) {
  install.packages("BiocManager", lib = site_lib, Ncpus = getOption("Ncpus", 2L))
}

# ----- Bioconductor package set -----
bioc_pkgs <- c(
  "BiocGenerics","Biobase","S4Vectors","IRanges",
  "Biostrings","GenomicRanges","GenomicFeatures",
  "AnnotationDbi","DESeq2","edgeR","limma",
  "SummarizedExperiment","SingleCellExperiment",
  "BiocParallel","tximport"
)

# Install Bioconductor first (will pull compatible versions)
BiocManager::install(bioc_pkgs, lib = site_lib, ask = FALSE,
                     update = TRUE, Ncpus = getOption("Ncpus", 2L))

# ----- Core CRAN set -----
cran_pkgs <- c(
  "tidyverse","ragg","data.table","readxl","janitor",
  "devtools","ggpubr","lubridate","forcats","remotes","lattice"
)

# Try current MASS/Matrix from CRAN first; if that fails, fall back to archived versions for R 4.3
safe_install <- function(pkgs, ...) {
  tryCatch(install.packages(pkgs, ...),
           error = function(e) { message("install.packages error: ", conditionMessage(e)) })
}

safe_install(c("MASS","Matrix"), lib = site_lib, Ncpus = getOption("Ncpus", 2L))
if (!requireNamespace("MASS", quietly = TRUE, lib.loc = site_lib))
  safe_install("https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-59.tar.gz",
               lib = site_lib, type = "source", Ncpus = getOption("Ncpus", 2L))
if (!requireNamespace("Matrix", quietly = TRUE, lib.loc = site_lib))
  safe_install("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz",
               lib = site_lib, type = "source", Ncpus = getOption("Ncpus", 2L))

# Now install the rest of CRAN
safe_install(cran_pkgs, lib = site_lib, Ncpus = getOption("Ncpus", 2L))

# ----- Summary -----
sel <- unique(c("BiocManager", bioc_pkgs, cran_pkgs, "MASS", "Matrix"))
ip <- installed.packages(lib.loc = site_lib)
keep <- rownames(ip) %in% sel
message("\nInstalled (site lib):")
print(ip[keep, c("Package","Version"), drop = FALSE])
RSCRIPT

# ---- Relock (chmod only) ----
echo "🔒 Relocking $R_LIBS_SITE..."
chmod -R a+rX "$R_LIBS_SITE"
chmod -R a-w  "$R_LIBS_SITE"

echo "✅ Installation complete. Shared library ready at $R_LIBS_SITE"

