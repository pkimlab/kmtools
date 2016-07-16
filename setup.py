from setuptools import setup, find_packages

setup(
    name='kmtools',
    version='0.0.7',
    author='kimlab.org',
    packages=find_packages(),
    package_data={
        'kmtools': [
            'cluster_tools/scripts/*.sh']},
)
