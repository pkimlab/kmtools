#!/bin/bash

set -ev

PKG_DIR="${RECIPE_DIR}/.."

python -m pytest \
    -c ./setup.cfg \
    --cov="${SP_DIR}/${PKG_NAME}/${SUBPKG_NAME}" \
    --benchmark-disable \
    "${PKG_DIR}/tests/${SUBPKG_NAME}"

mv .coverage "${PKG_DIR}/.coverage"
