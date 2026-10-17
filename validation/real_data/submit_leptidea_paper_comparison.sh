#!/usr/bin/env bash
set -euo pipefail
SCT_LP_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_LP_REPO"
SCT_LP_SOURCE=$(realpath "${1:-local_results/source_recovery_CRtNNj/leptidea}")
SCT_LP_TOOLS=$(realpath "${SCT_LEPTIDEA_TOOLS:-local_results/annotated_lep_rjv857/tools}")
SCT_LP_GFFREAD=$(realpath "${SCT_LEPTIDEA_GFFREAD:-local_results/rhynchospora_native_mNeWPw/tools/bin/gffread}")
for exe in python blastp makeblastdb MCScanX; do test -x "$SCT_LP_TOOLS/bin/$exe"; done
test -x "$SCT_LP_GFFREAD"
test -f "$SCT_LP_SOURCE/homology/COMPLETE"
grep -q -- '--skip-incomplete-variants' validation/benchmark/genespace_native_order.R
SCT_LP_RUN=$(mktemp -d "$SCT_LP_REPO/local_results/leptidea_paper_XXXXXX")
SCT_LP_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_LP_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
mkdir -p "$SCT_LP_RUN/code" logs
cp -a src validation "$SCT_LP_RUN/code/"
printf '%s\n' '{"datasets":[{"id":"leptidea_paper_assemblies"}]}' > "$SCT_LP_RUN/comparison_manifest.json"
printf '%s\n' "$SCT_LP_RUN" > local_results/latest_leptidea_paper_run.txt
export SCT_LP_REPO SCT_LP_SOURCE SCT_LP_TOOLS SCT_LP_GFFREAD SCT_LP_RUN SCT_LP_MAMBA SCT_LP_GS_ENV
SCT_LP_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time=00:30:00 --job-name=st_lep_protein --output=logs/st_lep_protein.%j.out --error=logs/st_lep_protein.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LP_RUN/code/src" PYTHONUNBUFFERED=1
"$SCT_LP_TOOLS/bin/python" "$SCT_LP_RUN/code/validation/real_data/leptidea_paper_comparison.py" prepare --run "$SCT_LP_RUN" --source "$SCT_LP_SOURCE" --gffread "$SCT_LP_GFFREAD"
SLURM
)
SCT_LP_PAIR=$(sbatch --parsable --partition=cpu --array=0-5 --cpus-per-task=4 --mem=8G --time=03:00:00 --dependency="afterok:${SCT_LP_SETUP%%;*}" --job-name=st_lep_paper_pair --output=logs/st_lep_paper_pair.%A_%a.out --error=logs/st_lep_paper_pair.%A_%a.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PATH="$SCT_LP_TOOLS/bin:$PATH" PYTHONPATH="$SCT_LP_RUN/code/src" PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
python "$SCT_LP_RUN/code/validation/real_data/leptidea_paper_comparison.py" pair --run "$SCT_LP_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_LP_FIG=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=00:30:00 --dependency="afterok:${SCT_LP_PAIR%%;*}" --job-name=st_lep_paper_fig --output=logs/st_lep_paper_fig.%j.out --error=logs/st_lep_paper_fig.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LP_RUN/code/src" PYTHONUNBUFFERED=1 MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$SCT_LP_TOOLS/bin/python" "$SCT_LP_RUN/code/validation/real_data/leptidea_paper_comparison.py" collect --run "$SCT_LP_RUN"
eval "$("$SCT_LP_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_LP_GS_ENV"
"$SCT_LP_TOOLS/bin/python" "$SCT_LP_RUN/code/validation/real_data/plot_comparison.py" --prepared "$SCT_LP_RUN/prepared/leptidea_paper_assemblies" --output "$SCT_LP_RUN/figures/leptidea_paper_assemblies" --seconds 300
"$SCT_LP_TOOLS/bin/python" "$SCT_LP_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_LP_RUN" --manifest "$SCT_LP_RUN/comparison_manifest.json" --latest-archive "$SCT_LP_REPO/local_results/SynTangle_Leptidea_paper_comparisons.zip"
SLURM
)
printf 'Protein preparation: %s\nSix pairwise jobs, no throttle: %s\nComparison and archive: %s\nResults: %s\nArchive when complete: %s/local_results/SynTangle_Leptidea_paper_comparisons.zip\n' "$SCT_LP_SETUP" "$SCT_LP_PAIR" "$SCT_LP_FIG" "$SCT_LP_RUN" "$SCT_LP_REPO"
