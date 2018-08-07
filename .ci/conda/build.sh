#!/bin/bash

set -ev

ln -s ${CC} $(dirname ${CC})/cc
ln -s ${GCC} $(dirname ${GCC})/gcc

${PYTHON} setup.py install --single-version-externally-managed --record=record.txt
