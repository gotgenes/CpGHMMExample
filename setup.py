#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path
import sys

from setuptools import setup

SCRIPTS = [os.path.sep.join(('scripts', script)) for script in
           os.listdir('scripts')]
REQUIREMENTS = []

setup(
    name='CpGHMMExample',
    version='0.1a1',
    author='Christopher D. Lasher',
    author_email='chris.lasher@gmail.com',
    install_requires=REQUIREMENTS,
    packages=['CpGHMMExample'],
    scripts=SCRIPTS,
    url='http://pypi.python.org/pypi/CpGHMMExample',
    license='MIT License',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    description=("Python bioinformatics tools for the Ovcharenko group"),
    long_description=open('README.rst').read(),
)

