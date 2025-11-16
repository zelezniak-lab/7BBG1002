#!/bin/bash -l
#SBATCH --job-name=remove-conda
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00

set -u

echo "=== Removing Conda from user account ==="

# --------------------------------------------------------
# 1) Deactivate any active Conda environment (if available)
# --------------------------------------------------------
if command -v conda >/dev/null 2>&1; then
  echo "Deactivating any active Conda environment..."
  # Deactivate repeatedly to escape nested activations
  for _ in 1 2 3; do conda deactivate >/dev/null 2>&1 || true; done
fi

# --------------------------------------------------------
# 2) Remove Conda paths from the current session PATH
# --------------------------------------------------------
if [[ ":$PATH:" == *"conda"* || ":$PATH:" == *"anaconda"* || ":$PATH:" == *"miniconda"* ]]; then
  echo "Removing Conda-related paths from current PATH..."
  NEWPATH="$(echo "$PATH" | tr ':' '\n' | grep -v -E '/(mini|ana)?conda[^/]*/bin' | grep -v -E '(mini|ana)?conda' | paste -sd ':' -)"
  export PATH="$NEWPATH"
fi

# --------------------------------------------------------
# 3) Candidates for removal
# --------------------------------------------------------
CONDA_DIRS=(
  "$HOME/miniconda3"
  "$HOME/anaconda3"
  "$HOME/.conda"
  "$HOME/.continuum"
  "$HOME/.condarc"
  "$HOME/.conda_environments.txt"
)

# --------------------------------------------------------
# 4) Remove Conda files/directories (if present)
# --------------------------------------------------------
for d in "${CONDA_DIRS[@]}"; do
  if [ -e "$d" ]; then
    echo "Removing $d ..."
    rm -rf -- "$d"
  fi
done

# --------------------------------------------------------
# 5) Clean shell initialisation files (remove full conda block)
# --------------------------------------------------------
timestamp="$(date +%Y%m%d-%H%M%S)"
SHELL_FILES=( "$HOME/.bashrc" "$HOME/.bash_profile" )

for f in "${SHELL_FILES[@]}"; do
  if [ -f "$f" ]; then
    echo "Cleaning Conda initialisation from $f ..."
    cp -a -- "$f" "${f}.bak.${timestamp}"

    # Remove the entire conda init block if present
    sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$f"

    # Remove any remaining stray lines referencing conda/anaconda/miniconda
    sed -i '/[Cc]onda/d' "$f"
    sed -i '/[Aa]naconda/d' "$f"
    sed -i '/[Mm]iniconda/d' "$f"

    # Trim possible blank-line spam left behind
    awk 'BEGIN{blank=0} {if ($0 ~ /^$/) {blank++} else {blank=0} if (blank<=1) print}' "$f" > "$f.tmp.$$" && mv "$f.tmp.$$" "$f"
  fi
done

# --------------------------------------------------------
# 6) Summary
# --------------------------------------------------------
echo "=== Done. Conda has been removed from your user environment. ==="
echo "Backups of edited shell files were saved with suffix: .bak.${timestamp}"
echo "Start a new shell or log out/in to ensure PATH is fully refreshed."

