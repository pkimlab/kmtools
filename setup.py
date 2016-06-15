from setuptools import setup, find_packages

setup(
    name='ascommon',
    version='0.0.6',
    author='Alexey Strokach',
    author_email='alex.strokach@utoronto.ca',
    packages=find_packages(),
    package_data={
        'ascommon': [
            'cluster_tools/scripts/*.sh']},
)
