#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path
import sys

from setuptools import setup
from setuptools.extension import Extension

SCRIPTS = [os.path.sep.join(('scripts', script)) for script in
           os.listdir('scripts')]
REQUIREMENTS = ['matplotlib', 'numpy']

setup(
        name='CpGHMMExample',
        version='0.1a1',
        author='Christopher D. Lasher',
        author_email='chris.lasher@gmail.com',
        install_requires=REQUIREMENTS,
        packages=['CpGHMMExample'],
        ext_modules=[Extension('CpGHMMExample.alglib', ['src/alglib.c'])],
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
        description=("An example of an HMM to detect CpG sites."),
        long_description=open('README.rst').read(),
)

