import os.path as op
import warnings

from setuptools import find_packages, setup


def _read_md_as_rst(file):
    """Read Markdown file and convert it to ReStructuredText."""
    from pypandoc import convert_file
    return convert_file(file, 'rst', format='md')


def _read_md_as_md(file):
    """Read Markdown file."""
    with open(op.join(op.dirname(__file__), file)) as ifh:
        return ifh.read()


def read_md(file):
    """Read MarkDown file and try to convert it to ReStructuredText if you can."""
    try:
        return _read_md_as_rst(file)
    except ImportError:
        warnings.warn("pypandoc module not found, could not convert Markdown to RST!")
        return _read_md_as_md(file)


setup(
    name='kmtools',
    version='0.0.26',
    author='kimlab.org',
    author_email='alex.strokach@utoronto.ca',
    url="https://gitlab.com/kimlab/kmtools",
    description="Bits of reusable code to make our lives easier.",
    long_description=read_md("README.md"),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license='MIT',
    packages=find_packages(),
    package_data={
        'kmtools.sequence_tools': ['support/*.csv'],
        'kmtools.structure_tools': ['data/*'],
    },
)
