from setuptools import setup, find_packages

setup(
    name='kmtools',
    version='0.0.13',
    author='kimlab.org',
    packages=['kmtools.' + x for x in find_packages('kmtools')],
    namespace_packages=['kmtools'],
    package_data={
        'kmtools': [
            'cluster_tools/scripts/*.sh',
            'sequence_tools/support/*.tsv']},
)
