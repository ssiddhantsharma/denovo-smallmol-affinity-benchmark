#!/usr/bin/env bash
# Score the GPU-only columns on the design box, unattended.
#
# WHY THIS IS NOT A GITHUB ACTION: hosted runners have no GPU, and a SELF-hosted
# runner on a PUBLIC repo lets any fork's pull request execute code on this
# machine -- GitHub says not to do it, and they are right. So this pulls instead
# of being pushed to: nothing inbound is exposed and no runner token exists.
#
# Install on the GPU box (NOT in CI):
#   crontab -e
#   30 10 * * 1  /path/to/repo/scripts/gpu_autoscore.sh >> ~/autoscore.log 2>&1
#
# Monday 10:30, an hour after the weekly watcher, so a PR merged that morning is
# picked up the same day. Requires: a conda/venv with the predictor deps, and
# push rights on the repo (an SSH deploy key is enough).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

# Overlap protection. A co-fold run can outlast the gap to the next cron, and
# two of them would fight over the same output dir and the same git index.
exec 9>"${TMPDIR:-/tmp}/gpu_autoscore.lock"
flock -n 9 || { echo "$(date -Is) another run holds the lock; exiting"; exit 0; }

echo "=== $(date -Is) autoscore ==="
git fetch --quiet origin
git checkout --quiet main
git pull --quiet --ff-only origin main   # ff-only: never silently merge local edits

n="$(python scripts/update.py --count)"
echo "systems needing scores: $n"
[ "$n" -eq 0 ] && { echo "nothing to do"; exit 0; }

# --cofold only. The paid Rowan half stays manual: this script runs unattended
# and spending money unattended on a box nobody is watching is a different
# decision from using a GPU that is already paid for.
python scripts/update.py --cofold --max-systems "${MAX_SYSTEMS:-20}"

python scripts/qc_dataset.py   # never push a set that fails its own checks

if git diff --quiet -- predictions/; then
  echo "scored nothing new"; exit 0
fi
git add predictions/
git commit -m "autoscore: GPU columns for $n system(s)"
git push origin main
echo "pushed; figures.yml will redraw"
