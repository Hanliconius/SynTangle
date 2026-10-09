#!/usr/bin/env bash
set -euo pipefail
if (( $# < 3 || $# > 4 )); then
    echo 'Usage: bash certification/submit.sh check|verify|export FIXTURE LAYOUT [CERTIFICATE_TO_VERIFY]' >&2
    exit 2
fi
SCT_PROOF_ROOT=$(cd "$(dirname "$0")/.." && pwd)
SCT_PROOF_ACTION=$1
SCT_PROOF_FIXTURE=$(realpath "$2")
SCT_PROOF_LAYOUT=$(realpath "$3")
case "$SCT_PROOF_ACTION" in check|verify|export) ;; *) exit 2 ;; esac
SCT_PROOF_RUN=$(mktemp -d "$SCT_PROOF_ROOT/proof_run.XXXXXX")
cp "$SCT_PROOF_ROOT/certification/proof.py" "$SCT_PROOF_RUN/proof.py"
cp "$SCT_PROOF_FIXTURE" "$SCT_PROOF_RUN/fixture.json"
cp "$SCT_PROOF_LAYOUT" "$SCT_PROOF_RUN/layout.json"
SCT_PROOF_OUTPUT="$SCT_PROOF_RUN/certificate.json"
if [[ "$SCT_PROOF_ACTION" == verify ]]; then
    [[ $# == 4 ]] || { echo 'verify requires certificate path' >&2; exit 2; }
    cp "$4" "$SCT_PROOF_OUTPUT"
elif [[ "$SCT_PROOF_ACTION" == export ]]; then
    SCT_PROOF_OUTPUT="$SCT_PROOF_RUN/model"
fi
export SCT_PROOF_RUN SCT_PROOF_ACTION SCT_PROOF_OUTPUT
sbatch --partition=cpu --job-name=st_exact_check --cpus-per-task=1 --mem=4G \
    --time=00:15:00 --output="$SCT_PROOF_RUN/slurm.%j.out" \
    --error="$SCT_PROOF_RUN/slurm.%j.err" --export=ALL <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
python "$SCT_PROOF_RUN/proof.py" "$SCT_PROOF_ACTION" \
    "$SCT_PROOF_RUN/fixture.json" "$SCT_PROOF_RUN/layout.json" "$SCT_PROOF_OUTPUT" \
    --seconds 600 --max-states 1000000
SLURM
printf 'Results and logs: %s\n' "$SCT_PROOF_RUN"
