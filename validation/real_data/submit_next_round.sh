#!/usr/bin/env bash
# Independent production snapshots and a separate certification checkout.
set -euo pipefail
SCT_NR_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_NR_REPO"
SCT_NR_PROOF=$(mktemp -d /gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_proof.XXXXXX)
git clone -q --single-branch --branch exact-certification-experiment https://github.com/Hanliconius/SynTangle.git "$SCT_NR_PROOF"
test -f "$SCT_NR_PROOF/certification/submit_flatworm_continuation.sh"
bash validation/real_data/submit_leptidea_followup.sh
bash validation/real_data/submit_codfish_carex.sh
bash "$SCT_NR_PROOF/certification/submit_flatworm_continuation.sh"
