#!/bin/bash

set -ev

PACKAGE_ROOT_DIR="${RECIPE_DIR}/../.."

python -m pytest \
    -c setup.cfg \
    --cov="${SP_DIR}/kmtools/${SUBPKG_NAME}" \
    --cov-config=setup.cfg \
    --benchmark-disable \
    --color=yes \
    "tests/${SUBPKG_NAME}"

# sed -i "s|${SP_DIR}||g" .coverage
mv .coverage "${PACKAGE_ROOT_DIR}/.coverage"
