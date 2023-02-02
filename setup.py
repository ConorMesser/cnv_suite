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
    version="0.1.1",
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
        'pandas>=1.1.3, <2',
        'numpy>=1.19.0, <1.23',
        'scipy>=1.5.3, <1.9',
        'matplotlib>=3.3.2, <3.6',
        'plotly>=5.6.0, <6',
        'intervaltree>=3.1.0, <4',
        'natsort>=8.1.0, <10',
        'pandarallel>=1.5.5, <1.8',
        'tqdm>=4.45',
        'kaleido>=0.2.1, <0.3',
        'nbformat>=4.2.0'
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