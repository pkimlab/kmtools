from setuptools import setup
from conda_build.metadata import MetaData

META = MetaData('devtools/conda-recipe')

setup(
    name=META.get_value('package/name'),
    version=META.get_value('package/version'),
    description=META.get_value('about/summary'),
    url=META.get_value('about/home'),
    author='Alexey Strokach',
    author_email='alex.strokach@utoronto.ca',
    packages=['ascommon'],
    # setup_requires=['pytest-runner'],
    # tests_require=['pytest'],
)
