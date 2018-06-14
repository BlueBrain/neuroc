#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

if sys.version_info < (2, 7):
    sys.exit("Sorry, Python < 2.7 is not supported")

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
        neuroc=neuroc.apps.__main__:cli
    ''',
    license="BBP-internal-confidential",
    install_requires=['numpy>=1.14.1',
                      'nose>=1.3.0',
                      'tqdm>=4.23.4',
                      'click==6.7',
                      'morphio>=0.9.9'],
    packages=find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
)
