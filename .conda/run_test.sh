#!/bin/bash

set -ev

PKG_DIR="${RECIPE_DIR}/.."

python -m pytest \
    -c "${PKG_DIR}/setup.cfg" \
    --cov="${SP_DIR}/${PKG_NAME}/${SUBPKG_NAME}" \
    --benchmark-disable \
    "${PKG_DIR}/${SUBPKG_NAME}"
