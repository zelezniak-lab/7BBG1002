#!/bin/bash -l
#SBATCH --job-name=ops-rstudio
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --signal=USR2
#SBATCH --cpus-per-task=1

# Must run as a batch job: the rstudio module only exists on the compute nodes.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "ERROR: submit this script with:  sbatch $0" 1>&2
  echo "       then read the connection details with:  cat slurm-<jobid>.out" 1>&2
  exit 1
fi

# ---- RStudio (AlmaLinux 10 / spack v1.2.1 stack) ----
# The rstudio module pulls in R plus the libraries needed to build tidyverse
# (freetype, libpng, libtiff, libjpeg-turbo, cairo, harfbuzz, fribidi,
# libxml2, curl, openssl, openblas, pkgconf, fontconfig).
module purge
module load rstudio/v2026.09.0_174 || { echo "ERROR: could not load the rstudio module" 1>&2; exit 1; }
module load cmake/3.31.11 libwebp/1.4.0

export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
export PKG_CONFIG_ALLOW_SYSTEM_LIBS=1

# ---- Runtime linker paths via pkg-config ----
for pc in freetype2 libpng libtiff-4 libjpeg libwebp libwebpmux; do
  d=$(pkg-config --variable=libdir "$pc" 2>/dev/null || true)
  [[ -n "$d" ]] && export LD_LIBRARY_PATH="$d${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
done

# ---- Shared + per-student R libraries (no writes to HOME) ----
# Libraries are tied to the R minor version, so derive it from R itself.
export R_VERSION_SHORT=$(Rscript -e 'cat(R.version$major, strsplit(R.version$minor, ".", fixed=TRUE)[[1]][1], sep=".")')
export SHARED_RLIB="/scratch/grp/msc_appbio/shared/Rlib/${R_VERSION_SHORT}"     # read-only for students
export CLASS_USER_ROOT="/scratch/grp/msc_appbio/userlibs/${R_VERSION_SHORT}"    # per-user, writable
export R_LIBS_USER="${CLASS_USER_ROOT}/${USER}"                                 # per-user lib (NOT HOME)
if [[ -d "$SHARED_RLIB" ]]; then
  export R_LIBS_SITE="$SHARED_RLIB"                                             # site library for R
  export R_LIBS="$SHARED_RLIB"                                                  # ensure it's on .libPaths()
else
  echo "WARNING: shared R library $SHARED_RLIB not found; only your personal library will be used." 1>&2
  unset R_LIBS R_LIBS_SITE
fi
umask 002
mkdir -p "$R_LIBS_USER" || { echo "ERROR: cannot create $R_LIBS_USER" 1>&2; exit 1; }
# ---------------------------------------------------------------------

# ---- Per-job working directory (TMPDIR is private to the job) ----
JOBTMP="${TMPDIR:-/tmp}/rstudio-${SLURM_JOB_ID:-$$}"
mkdir -p "$JOBTMP"

# ---- R options (mirror) ----
export R_PROFILE_USER="$JOBTMP/Rprofile.site"
cat > "$R_PROFILE_USER" <<'RPROF'
options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk"))
# Ensure our library order: user (writable) -> shared site -> system
user_lib <- Sys.getenv("R_LIBS_USER")
site_lib <- Sys.getenv("R_LIBS_SITE")
.libPaths(c(if (nzchar(user_lib)) user_lib, if (nzchar(site_lib)) site_lib, .Library.site, .Library))
RPROF
# ---------------------------------------------------------------------

# ---- RStudio server launch ----
export PASSWORD=$(openssl rand -base64 15)
readonly PORT=$(python3 -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
readonly NODE=$(hostname -s)
cat 1>&2 <<END
1. SSH tunnel from your workstation:

   Linux/Mac:
   ssh -NL 8787:${NODE}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   Windows:
   ssh -m hmac-sha2-512 -NL 8787:${NODE}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   Then open http://localhost:8787

2. Login to RStudio:
   user: ${USER}
   pass: ${PASSWORD}

To end: close the RStudio session and run:  scancel -f ${SLURM_JOB_ID}

To SSH to the compute node directly:  ssh ${NODE}
END

DBCONF="$JOBTMP/database.conf"
echo "directory=$JOBTMP/var-rstudio-server" > "$DBCONF"

# ---- rsession wrapper ----
# rserver starts each R session through PAM with a clean environment, so the
# library paths, R profile and linker paths set above never reach R. This
# wrapper restores them before starting the real rsession.
RSESSION_WRAPPER="$JOBTMP/rsession.sh"
{
  echo '#!/bin/bash'
  for v in R_LIBS_USER R_LIBS_SITE R_LIBS R_PROFILE_USER LD_LIBRARY_PATH PATH \
           PKG_CONFIG_PATH PKG_CONFIG_ALLOW_SYSTEM_CFLAGS PKG_CONFIG_ALLOW_SYSTEM_LIBS; do
    [[ -n "${!v:-}" ]] && printf 'export %s=%q\n' "$v" "${!v}"
  done
  # Make `module` (and so `module load anaconda3/...`) work in RStudio's
  # Terminal and in notebook bash chunks.
  echo 'source /etc/profile.d/99-environment.sh'
  echo 'source /etc/profile.d/01-lmod.sh'
  echo 'export -f module ml 2>/dev/null'
  printf 'exec %q "$@"\n' "$(command -v rsession)"
} > "$RSESSION_WRAPPER"
chmod 700 "$RSESSION_WRAPPER"
# ---------------------------------------------------------------------

rserver --rsession-path="$RSESSION_WRAPPER" --server-user "${USER}" --www-port "${PORT}" --server-data-dir "$JOBTMP/data-rstudio-server" \
  --secure-cookie-key-file "$JOBTMP/data-rstudio-server/secure-cookie-key" \
  --database-config-file="$DBCONF" --auth-none=0 \
  --auth-pam-helper-path=pam-env-helper

printf 'RStudio Server exited\n' 1>&2
