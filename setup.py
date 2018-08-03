import os.path as op

from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup


def read_md(file):
    """Read Markdown file."""
    with open(op.join(op.dirname(__file__), file)) as ifh:
        return ifh.read()


EXTENSIONS = [Extension("*", ["kmtools/sequence_tools/*.pyx"])]


setup(
    name="kmtools",
    version="0.0.27",
    author="kimlab.org",
    author_email="alex.strokach@utoronto.ca",
    url="https://gitlab.com/kimlab/kmtools",
    description="Bits of reusable code to make our lives easier.",
    long_description=read_md("README.md"),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
    packages=find_packages(),
    package_data={
        "kmtools.sequence_tools": ["support/*.csv"],
        "kmtools.structure_tools": ["data/*"],
    },
    ext_modules=cythonize(EXTENSIONS),
)
