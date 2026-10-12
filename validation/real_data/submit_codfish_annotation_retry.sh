#!/usr/bin/env bash
set -euo pipefail
SCT_CD_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_CD_REPO"
SCT_CD_RUN=$(realpath "${1:-local_results/codfish_carex_kdt7Zz/codfish}")
SCT_CD_PAIR=${2:-138436841}
SCT_CD_FIG=${3:-138436842}
SCT_CD_TOOLS=$(realpath local_results/annotated_lep_rjv857/tools)
SCT_CD_GFFREAD=$(realpath local_results/rhynchospora_native_mNeWPw/tools/bin/gffread)
test -x "$SCT_CD_TOOLS/bin/python"
test -x "$SCT_CD_GFFREAD"
for name in Melanogrammus_aeglefinus Boreogadus_saida Arctogadus_glacialis Lota_lota Merluccius_merluccius; do
  test -f "$SCT_CD_RUN/sources/$name/COMPLETE"
done
SCT_CD_PAIR_STATE=$(scontrol show job -o "$SCT_CD_PAIR")
SCT_CD_FIG_STATE=$(scontrol show job -o "$SCT_CD_FIG")
[[ "$SCT_CD_PAIR_STATE" == *'JobName=st_codfish_pair '* && "$SCT_CD_PAIR_STATE" == *'JobState=PENDING '* ]]
[[ "$SCT_CD_FIG_STATE" == *'JobName=st_codfish_fig '* && "$SCT_CD_FIG_STATE" == *'JobState=PENDING '* ]]
SCT_CD_REPAIR=$(mktemp -d "$SCT_CD_RUN/annotation_retry_XXXXXX")
mkdir -p "$SCT_CD_REPAIR/code" logs
cp -a src validation "$SCT_CD_REPAIR/code/"
export SCT_CD_RUN SCT_CD_REPAIR SCT_CD_TOOLS SCT_CD_GFFREAD
SCT_CD_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time=00:30:00 --job-name=st_codfish_retry --output=logs/st_codfish_retry.%j.out --error=logs/st_codfish_retry.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_CD_REPAIR/code/src:$SCT_CD_REPAIR/code/validation/real_data" PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
"$SCT_CD_TOOLS/bin/python" -m unittest test_paper_genome_reanalysis > "$SCT_CD_REPAIR/tests.log" 2>&1
"$SCT_CD_TOOLS/bin/python" "$SCT_CD_REPAIR/code/validation/real_data/paper_genome_reanalysis.py" prepare --run "$SCT_CD_RUN" --manifest "$SCT_CD_REPAIR/code/validation/real_data/codfish_carex_manifest.json" --case 0 --index 0 --gffread "$SCT_CD_GFFREAD"
SLURM
)
printf 'Atlantic cod repair: %s\nAudit: %s/sources/Gadus_morhua_NEAC/sequence_mapping_audit.json\n' "$SCT_CD_JOB" "$SCT_CD_RUN"
scontrol update JobId="$SCT_CD_PAIR" Dependency="afterok:${SCT_CD_JOB%%;*}"
scontrol update JobId="$SCT_CD_FIG" Dependency="afterok:$SCT_CD_PAIR"
printf 'Reconnected existing pair array %s and comparison %s; other preparations reused.\n' "$SCT_CD_PAIR" "$SCT_CD_FIG"
