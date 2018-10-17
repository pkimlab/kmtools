# kmtools

[![anaconda](https://img.shields.io/conda/dn/kimlab/kmtools.svg)](https://anaconda.org/kimlab/kmtools/)
[![docs](https://img.shields.io/badge/docs-v0.0.35-blue.svg?version=latest)](https://kimlab.gitlab.io/kmtools/v0.0.35/)
[![build status](https://gitlab.com/kimlab/kmtools/badges/v0.0.35/build.svg)](https://gitlab.com/kimlab/kmtools/commits/v0.0.35/)
[![coverage report](https://gitlab.com/kimlab/kmtools/badges/v0.0.35/coverage.svg)](https://kimlab.gitlab.io/kmtools/v0.0.35/htmlcov/)

Bits of reusable code to make our lives easier.

Follows the *Whole Sort of General Mish Mash* design principle.

## Contents

- [kmtools](#kmtools)
    - [Contents](#contents)
    - [Tools](#tools)
        - [Structure tools](#structure-tools)
            - [To do](#to-do)
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
