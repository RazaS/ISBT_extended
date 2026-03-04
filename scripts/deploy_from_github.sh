#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK_FILE="$REPO_DIR/.automation.lock"

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  echo "[autodeploy] lock busy; skipping"
  exit 0
fi

cd "$REPO_DIR"

git fetch origin main

LOCAL_SHA="$(git rev-parse HEAD)"
REMOTE_SHA="$(git rev-parse origin/main)"

if [[ "$LOCAL_SHA" == "$REMOTE_SHA" ]]; then
  echo "[autodeploy] already up to date"
  exit 0
fi

if [[ -n "$(git status --porcelain)" ]]; then
  echo "[autodeploy] working tree is dirty; skipping deploy"
  exit 0
fi

if ! git merge-base --is-ancestor HEAD origin/main; then
  echo "[autodeploy] local branch has commits not in origin; skipping deploy"
  exit 0
fi

git pull --ff-only origin main

if [[ ! -x "$REPO_DIR/.venv/bin/python" ]]; then
  python3 -m venv "$REPO_DIR/.venv"
fi

"$REPO_DIR/.venv/bin/pip" install -q -r requirements.txt

systemctl restart isbt-extended.service

echo "[autodeploy] deployed $REMOTE_SHA"
