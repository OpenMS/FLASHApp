#!/usr/bin/env bash
# FLASHApp's step of the desktop build (deploy.desktop_app.stage_script in
# app.yaml), run by build-desktop-app.yaml once requirements.txt is installed:
#   bash scripts/desktop-stage.sh <python of the payload> <payload dir>
set -euo pipefail
py=$1

"$py" -c 'import openms_insight; print("openms_insight", openms_insight.__version__)'

case "$(uname -s)" in
  Darwin) ;;
  *)
    # polars' default wheel needs AVX2 and more; the LTS build runs on any
    # x86_64 CPU. Both ship the polars/ package, so remove the default one
    # first: installing over it mixed files of two versions (ImportError in
    # polars.io).
    "$py" -m pip uninstall -y --disable-pip-version-check polars \
      $("$py" -m pip list --format=freeze --disable-pip-version-check | grep -o '^polars-runtime[^=]*' || true)
    "$py" -m pip install --no-warn-script-location --disable-pip-version-check 'polars-lts-cpu >= 1.0, < 2'
    "$py" -c 'import polars, polars.io.delta; print("polars", polars.__version__)'
    ;;
esac
