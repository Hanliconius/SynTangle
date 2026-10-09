#!/usr/bin/env bash
set -euo pipefail
if (( $# < 3 || $# > 4 )); then
    echo 'Usage: bash certification/submit_structural.sh check|verify FIXTURE LAYOUT [CERTIFICATE_TO_VERIFY]' >&2
    exit 2
fi
SCT_STRUCT_ROOT=$(cd "$(dirname "$0")/.." && pwd)
SCT_STRUCT_ACTION=$1
case "$SCT_STRUCT_ACTION" in check|verify) ;; *) exit 2 ;; esac
SCT_STRUCT_RUN=$(mktemp -d "$SCT_STRUCT_ROOT/structural_run.XXXXXX")
cp "$SCT_STRUCT_ROOT/certification/proof.py" "$SCT_STRUCT_RUN/proof.py"
cp "$SCT_STRUCT_ROOT/certification/structural.py" "$SCT_STRUCT_RUN/structural.py"
cp "$2" "$SCT_STRUCT_RUN/fixture.json"
cp "$3" "$SCT_STRUCT_RUN/layout.json"
if [[ "$SCT_STRUCT_ACTION" == verify ]]; then
    [[ $# == 4 ]] || { echo 'verify requires certificate path' >&2; exit 2; }
    cp "$4" "$SCT_STRUCT_RUN/certificate.json"
fi
export SCT_STRUCT_RUN SCT_STRUCT_ACTION
sbatch --partition=cpu --job-name=st_structural --cpus-per-task=1 --mem=4G \
    --time=00:15:00 --output="$SCT_STRUCT_RUN/slurm.%j.out" \
    --error="$SCT_STRUCT_RUN/slurm.%j.err" --export=ALL <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
python "$SCT_STRUCT_RUN/structural.py" "$SCT_STRUCT_ACTION" \
    "$SCT_STRUCT_RUN/fixture.json" "$SCT_STRUCT_RUN/layout.json" "$SCT_STRUCT_RUN/certificate.json"
SLURM
printf 'Results and logs: %s\n' "$SCT_STRUCT_RUN"
