#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK_FILE="$REPO_DIR/.automation.lock"

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  echo "[isbt-refresh] lock busy; skipping"
  exit 0
fi

cd "$REPO_DIR"

if [[ ! -x "$REPO_DIR/.venv/bin/python" ]]; then
  python3 -m venv "$REPO_DIR/.venv"
fi

# Keep local checkout aligned with remote before generating new data.
if [[ -z "$(git status --porcelain)" ]]; then
  git pull --rebase origin main
else
  echo "[isbt-refresh] working tree is dirty; skipping pull"
fi

STATUS="$($REPO_DIR/.venv/bin/python "$REPO_DIR/scripts/refresh_isbt_cache.py")"
if [[ "$STATUS" != "changed" ]]; then
  echo "[isbt-refresh] no upstream ISBT data changes"
  exit 0
fi

git add isbt_cache.json
if git diff --cached --quiet; then
  echo "[isbt-refresh] no git diff after refresh"
  exit 0
fi

git config user.name "${GIT_BOT_NAME:-ISBT Bot}"
git config user.email "${GIT_BOT_EMAIL:-isbt-bot@bloodapps.com}"

STAMP="$(date -u +"%Y-%m-%d %H:%M UTC")"
git commit -m "Refresh ISBT cache (${STAMP})"

git push origin main

# Restart app so freshly committed cache is what the app serves.
systemctl restart isbt-extended.service

echo "[isbt-refresh] cache refreshed, committed, pushed, and service restarted"
