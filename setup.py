#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

VERSION = imp.load_source("", "neuroc/version.py").__version__

setup(
    name="neuroc",
    author="BlueBrain NSE",
    author_email="bbp-ou-nse@groupes.epfl.ch",
    version=VERSION,
    description='NeuroC: a collection of tools for morphology cloning applications',
    url="https://bbpteam.epfl.ch/project/issues/projects/NSETM/issues",
    download_url="ssh://bbpcode.epfl.ch/nse/morph-tool",
    entry_points='''
        [console_scripts]
        neuroc=neuroc.cli:cli
    ''',
    license="BBP-internal-confidential",
    python_requires='>=3.6',
    install_requires=['attrs>=19.1.0',
                      'numpy>=1.15.1',
                      'tqdm>=4.23.4',
                      'click>=6.7',
                      'pandas>=1.0.3',
                      'morphio>=3.0,<4.0',
                      'morph-tool>=2.9.0,<3.0',
                      'neurom>=3.0,<4.0',
                      'scikit-learn>=0.21.3',
    ],
    extras_require={
        'plotly': [
            'neurom[plotly]>=3.0,<4.0',
            'dash>=1.1.1',  # The core dash backend
            'dash-html-components>=1.0.0',  # HTML components
            'dash-core-components>=1.1.1',  # Supercharged components
            'dash-table>=4.1.0',  # Interactive DataTable component (new!)
        ]
    },
    packages=find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
)
