# Deployment Automation

This repo includes two automation scripts used on the VPS:

- `scripts/deploy_from_github.sh`
  - Pulls latest `origin/main` changes and restarts `isbt-extended.service`.
- `scripts/isbt_refresh_and_push.sh`
  - Forces a fresh ISBT API pull, updates `isbt_cache.json` only when substantive data changes,
    commits and pushes to `main`, then restarts `isbt-extended.service`.

Expected server path: `/opt/isbt_extended`

Both scripts use a shared lock file `.automation.lock` to prevent overlap.
