# KMTools

[![anaconda](https://anaconda.org/kimlab/kmtools/badges/version.svg?style=flat-square)](https://anaconda.org/kimlab/kmtools)
[![docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square&?version=latest)](https://kimlab.gitlab.org/kmtools)
[![build status](https://gitlab.com/kimlab/kmtools/badges/master/build.svg)](https://gitlab.com/kimlab/kmtools/commits/master)
[![coverage report](https://gitlab.com/kimlab/kmtools/badges/master/coverage.svg)](https://gitlab.com/kimlab/kmtools/commits/master)

Bits of reusable code to make our lives easier.

Follows the *Whole Sort of General Mish Mash* design principle.

## Contents

- [Tools](#tools)
  - [DB tools](#db-tools)
  - [DF tools](#df-tools)
  - [PY tools](#py-tools)
  - [Sequence tools](#sequence-tools)
  - [Structure tools](#structure-tools)
  - [System tools](#system-tools)
- [Contributing](#contributing)

## Tools

### Structure tools

Using [kmbio](https://github.com/kimlaborg/kmbio) instead of [biopython](https://github.com/biopython/biopython) leads to substantially better performance (> 2x with lots of room for improvement).

#### To do

- [ ] Cythonize more bottlenecks to improve performance.
- [ ] Add Python / Cython code for generating biological assemblies.

## Contributing

- Make sure all tests pass before merging into master.
- Follow the PEP8 / PyFlake / Flake8 / etc. guidelines.
- Add tests for new code.
- Try to document things.
- Break any / all of the above if you have a good reason.
