#!/bin/bash
# SessionStart hook: install FLASHApp deps + the OpenMS-Insight sibling repo so
# tests and the migration work in Claude Code on the web.
#
# Layout assumed (Claude Code on the web clones siblings under the same parent):
#   <parent>/FLASHApp                      (this repo, $CLAUDE_PROJECT_DIR)
#   <parent>/OpenMS-Insight                (visualization library dependency)
#
# The OpenMS-Insight wheel force-includes its built js-component/dist, so the
# Vue bundle must be built BEFORE the editable install. We handle that ordering
# and degrade gracefully if the sibling repo isn't present.
set -euo pipefail

# Only run in the remote (web) environment; local setups are managed by the user.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"
PARENT_DIR="$(dirname "$PROJECT_DIR")"
OI_DIR="$PARENT_DIR/OpenMS-Insight"

echo "[session-start] FLASHApp setup starting"

# --- 1. Python dependencies for FLASHApp -----------------------------------
# requirements.txt pins the app deps (streamlit, polars, pyopenms, scipy, ...).
# It also lists `openms-insight @ git+...`, but building that from the git URL
# fails (the wheel force-includes a js bundle that isn't built in a fresh
# clone) AND aborts the whole install. So strip that line here and install the
# local sibling separately in step 2.
REQ_TMP="$(mktemp)"
grep -ivE '^openms-insight[[:space:]]*@' "$PROJECT_DIR/requirements.txt" > "$REQ_TMP" || cp "$PROJECT_DIR/requirements.txt" "$REQ_TMP"
python3 -m pip install --user --quiet -r "$REQ_TMP" || {
  echo "[session-start] WARNING: requirements.txt install hit an error; continuing"
}
rm -f "$REQ_TMP"

# pytest for the test suite (not in the app requirements).
python3 -m pip install --user --quiet pytest pytest-cov

# --- 2. OpenMS-Insight sibling (build Vue bundle, then editable install) ----
if [ -d "$OI_DIR" ]; then
  echo "[session-start] Building OpenMS-Insight Vue bundle"
  # npm install is cache-friendly and idempotent (prefer over npm ci so the
  # cached container state is reused across sessions).
  ( cd "$OI_DIR/js-component" && npm install --no-audit --no-fund --silent && npm run build ) || {
    echo "[session-start] WARNING: Vue bundle build failed; OI install may be degraded"
  }

  # The wheel force-includes openms_insight/js-component/dist. Vite builds to
  # js-component/dist (repo root), so mirror it into the package path the build
  # backend expects, making the editable install succeed.
  if [ -d "$OI_DIR/js-component/dist" ]; then
    mkdir -p "$OI_DIR/openms_insight/js-component"
    rm -rf "$OI_DIR/openms_insight/js-component/dist"
    cp -r "$OI_DIR/js-component/dist" "$OI_DIR/openms_insight/js-component/dist"
  fi

  echo "[session-start] Installing OpenMS-Insight (editable)"
  python3 -m pip install --user --quiet -e "$OI_DIR" || {
    echo "[session-start] WARNING: editable OI install failed; falling back to PYTHONPATH"
    # Fallback: make it importable directly from source.
    echo "export PYTHONPATH=\"$OI_DIR:\${PYTHONPATH:-}\"" >> "${CLAUDE_ENV_FILE:-/dev/null}"
  }
else
  echo "[session-start] OpenMS-Insight not found at $OI_DIR (skipping; clone it as a sibling for the migration)"
fi

# --- 3. Persist environment for the session --------------------------------
# FLASHApp imports as top-level `src.*`; ensure the repo root is importable.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo "export PYTHONPATH=\"$PROJECT_DIR:\${PYTHONPATH:-}\"" >> "$CLAUDE_ENV_FILE"
  # User-site installs land here; make sure they're on PATH for pytest etc.
  echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$CLAUDE_ENV_FILE"
fi

echo "[session-start] FLASHApp setup complete"
