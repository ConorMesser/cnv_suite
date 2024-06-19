from setuptools import setup, find_packages
import sys

ver_info = sys.version_info
if ver_info < (3,7,0):
    raise RuntimeError("CNV_Suite requires at least python 3.7")

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='cnv_suite',
    version="0.2.0",
    packages=['cnv_suite.compare', 'cnv_suite.simulate', 'cnv_suite.utils', 'cnv_suite.visualize', 'cnv_suite'],
    entry_points={
        'console_scripts': [
            'simulate = cnv_suite.simulate.cnv_profile:main',
            'visualize = cnv_suite.visualize.plot_cnv_profile:main',
            'compare = cnv_suite.compare.__main__:main'
        ]
    },
    description='Copy Number tools for visualization, simulation, and comparison.',
    author='Conor Messer',
    author_email='cmesser@broadinstitute.org',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    install_requires=[
        'pandas',
        'numpy>=1.16',
        'scipy',
        'matplotlib>=3.1',
        'plotly>=5.12, <6',  # unknown if there will be breaking changes in a future major version of plotly
        'intervaltree',
        'natsort',
        'pandarallel',
        'tqdm',
        'kaleido',
        'nbformat',
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License"
    ],
    license="BSD3"
)
