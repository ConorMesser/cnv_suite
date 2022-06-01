from setuptools import setup
import sys

ver_info = sys.version_info
if ver_info < (3,7,0):
    raise RuntimeError("CNV_Suite requires at least python 3.7")

setup(
    name='CNV_Suite',
    version="0.1",
    packages=[
        'cnv_suite',
    ],
    # entry_points={
    #     'console_scripts': [
    #         'hapaseg = hapaseg.__main__:main',
    #     ]
    # },
    description='Copy Number tools for visualization, simulation, and comparison.',
    author='Conor Messer',
    author_email='cmesser@broadinstitute.org',
    # long_description = long_description,
    # long_description_content_type = 'text/markdown',
    install_requires=[
        'pandas>=1.1.3',
        'numpy>=1.19.0',
        'scipy>=1.5.3',
        'matplotlib>=3.3.2',
        'plotly>=5.6.0',
        'intervaltree>=3.1.0',
        'natsort>=8.1.0',
        'pandarallel>=1.5.5'
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