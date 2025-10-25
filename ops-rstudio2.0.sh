#!/bin/bash -l
#SBATCH --job-name=ops-rstudio
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --signal=USR2
#SBATCH --cpus-per-task=1

# ---- RStudio + stuff needed to install tidyverse ----
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

# ---- Runtime linker paths via pkg-config ----
for pc in freetype2 libpng libtiff-4 libjpeg libwebp libwebpmux; do
  d=$(pkg-config --variable=libdir "$pc" 2>/dev/null || true)
  [[ -n "$d" ]] && export LD_LIBRARY_PATH="$d:$LD_LIBRARY_PATH"
done

# ---- Shared + per-student R libraries (no writes to HOME) ----
export R_VERSION_SHORT=4.3
export SHARED_RLIB="/scratch/grp/msc_appbio/shared/Rlib/${R_VERSION_SHORT}"     # read-only for students
export CLASS_USER_ROOT="/scratch/grp/msc_appbio/userlibs/${R_VERSION_SHORT}"    # per-user, writable
export R_LIBS_SITE="$SHARED_RLIB"                                               # site library for R
export R_LIBS="$SHARED_RLIB${R_LIBS:+:$R_LIBS}"                                 # ensure it's on .libPaths()
export R_LIBS_USER="${CLASS_USER_ROOT}/${USER}"                                 # per-user lib (NOT HOME)
umask 002
mkdir -p "$R_LIBS_USER"
# ---------------------------------------------------------------------

# ---- R options (mirror) ----
export R_DEFAULT_PACKAGES=
export R_PROFILE_USER="${R_PROFILE_USER:-$TMPDIR/Rprofile.site}"
cat > "$R_PROFILE_USER" <<'RPROF'
options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk"))
# Ensure our library order: user (writable) -> shared site -> system
user_lib <- Sys.getenv("R_LIBS_USER")
site_lib <- Sys.getenv("R_LIBS")
if (nzchar(user_lib) && nzchar(site_lib)) {
  .libPaths(c(user_lib, site_lib, .Library.site, .Library))
}
RPROF
# ---------------------------------------------------------------------

# ---- RStudio server launch ----
export PASSWORD=$(openssl rand -base64 15)
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation:

   Linux/Mac:
   ssh -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   Windows:
   ssh -m hmac-sha2-512 -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   Then open http://localhost:8787

2. Login to RStudio:
   user: ${USER}
   pass: ${PASSWORD}

To end: close the RStudio session and run:  scancel -f ${SLURM_JOB_ID}

To SSH to the compute node directly:  ssh ${HOSTNAME}
END

DBCONF=$TMPDIR/database.conf
if [ ! -e "$DBCONF" ]; then
  printf "\nNOTE: creating $DBCONF database config file.\n\n"
  echo "directory=$TMPDIR/var-rstudio-server" > "$DBCONF"
fi

rserver --server-user "${USER}" --www-port "${PORT}" --server-data-dir "$TMPDIR/data-rstudio-server" \
  --secure-cookie-key-file "$TMPDIR/data-rstudio-server/secure-cookie-key" \
  --database-config-file="$DBCONF" --auth-none=0 \
  --auth-pam-helper-path=pam-env-helper

printf 'RStudio Server exited\n' 1>&2

