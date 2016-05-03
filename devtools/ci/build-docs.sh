#!/bin/bash

set -ev

if [[ -z "${GITHUB_TOKEN}" ]] ; then
    echo "GitHub API key needs to be set to update docs."
    exit -1
fi
cd "${TRAVIS_BUILD_DIR}/docs"

GITHUB_USERNAME=${TRAVIS_REPO_SLUG%/*}

# Install sphinx requirements
pip install -U pip
# Work around missing dependency in conda.
pip install auxlib
pip install -r requirements.txt

# Build the documentation
git config user.name "Travis CI"
mkdir -p "${PKG_NAME}.github.io"
cd "${PKG_NAME}.github.io"
git init
git checkout -b gh-pages
git remote add origin https://github.com/${GITHUB_USERNAME}/${PKG_NAME}.git
sphinx-build .. .
touch .nojekyll
git add .nojekyll
echo '.*' >> .gitignore
git add .
git commit --all -m "Updated docs to commit ${TRAVIS_COMMIT}."
git push -f -q "https://${GITHUB_TOKEN}@github.com/${GITHUB_USERNAME}/${PKG_NAME}.git" gh-pages \
    &> /dev/null

# Cleanup
cd ..
rm -rf "${PKG_NAME}.github.io"
