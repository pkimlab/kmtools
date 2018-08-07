# Installation

## Stable release

To install `kmtools`, run this command in your terminal:

```bash
conda install kmtools -c kimlab
```

This is the preferred method to install `kmtools`, as it will always install the most recent stable release.

If you don't have [conda] installed, this [Python installation guide] can guide
you through the process.

[conda]: https://conda.io
[Python installation guide]: https://conda.io/docs/user-guide/install/index.html

## From sources

The sources for `kmtools` can be downloaded from the [GitLab repo].

You can either clone the public repository:

```bash
git clone git://gitlab.com/kimlab/kmtools
```

Or download the [tarball]:

```bash
curl -OL https://gitlab.com/kimlab/kmtools/repository/master/archive.tar
```

Once you have a copy of the source, you can install it with:

```bash
python setup.py install
```

[GitLab repo]: https://gitlab.com/kimlab/kmtools
[tarball]: https://gitlab.com/kimlab/kmtools/repository/master/archive.tar
