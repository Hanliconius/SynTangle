#!/usr/bin/env bash
set -euo pipefail
SCT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_ROOT"
SCT_LIFT="$SCT_ROOT/local_results/leptidea_liftoff_0aWZxo"
SCT_TE="$SCT_ROOT/local_results/rhynchospora_TE_iBUOIo"
test -x "$SCT_LIFT/tools/bin/python"
SCT_RUN=$(mktemp -d "$SCT_ROOT/local_results/source_recovery_XXXXXX")
mkdir -p "$SCT_RUN/code" logs
cp -a validation/real_data "$SCT_RUN/code/"
export SCT_RUN SCT_LIFT SCT_TE
SCT_JOB=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=4G --time=00:25:00 --job-name=st_source_recovery --output="$SCT_ROOT/logs/st_source_recovery.%j.out" --error="$SCT_ROOT/logs/st_source_recovery.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONUNBUFFERED=1
cd "$SCT_RUN/code/real_data"
"$SCT_LIFT/tools/bin/python" -m unittest test_leptidea_recovery > "$SCT_RUN/leptidea_tests.log" 2>&1
"$SCT_LIFT/tools/bin/python" recover_leptidea_chromosomes.py --prior "$SCT_LIFT" --output "$SCT_RUN/leptidea"
SLURM
)
SCT_ZIP_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=01:00:00 --job-name=st_repeat_inventory --output="$SCT_ROOT/logs/st_repeat_inventory.%j.out" --error="$SCT_ROOT/logs/st_repeat_inventory.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONUNBUFFERED=1
cd "$SCT_RUN/code/real_data"
"$SCT_LIFT/tools/bin/python" inspect_repeat_archive.py --prior "$SCT_TE" --output "$SCT_RUN/repeats"
SLURM
)
printf 'Leptidea recovery: %s\nRepeat archive inspection: %s\nResults: %s\n' "$SCT_JOB" "$SCT_ZIP_JOB" "$SCT_RUN"
