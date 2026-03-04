#!/usr/bin/env python3
"""Refresh ISBT cache JSON by forcing a fresh pull from ISBT APIs.

Prints exactly one token to stdout:
- changed
- unchanged
"""

import copy
import importlib.util
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
APP_PATH = REPO_ROOT / "geno_mock_python" / "app.py"
CACHE_PATH = REPO_ROOT / "isbt_cache.json"


def load_app_module():
    spec = importlib.util.spec_from_file_location("isbt_app", APP_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load app module from {APP_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_cached_dataset(path: Path):
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def canonicalize(dataset):
    if not isinstance(dataset, dict):
        return dataset
    data = copy.deepcopy(dataset)
    meta = data.get("metadata")
    if isinstance(meta, dict):
        meta.pop("updated_at", None)
        meta.pop("updated_epoch", None)
    return data


def stable_json_bytes(payload):
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def write_dataset(path: Path, payload):
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False)
    tmp_path.replace(path)


def main():
    app = load_app_module()

    # Force a fresh pull (does not use the 30-day cache gate in load_isbt_dataset).
    fresh = app.build_isbt_dataset()
    cached = load_cached_dataset(CACHE_PATH)

    if cached is not None:
        if stable_json_bytes(canonicalize(cached)) == stable_json_bytes(canonicalize(fresh)):
            print("unchanged")
            return 0

    write_dataset(CACHE_PATH, fresh)
    print("changed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"refresh_isbt_cache.py failed: {exc}", file=sys.stderr)
        raise
