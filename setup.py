#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Setup for ejplugins."""

import io
from importlib import import_module
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()
with open('test_requirements.txt') as f:
    test_requirements = f.read().splitlines()

with io.open('README.md') as readme:
    setup(
        name='ejplugins',
        version=import_module('ejplugins').__version__,
        description='parser plugins for jsonextended',
        long_description=readme.read(),
        long_description_content_type='text/markdown',
        install_requires=requirements,
        tests_require=test_requirements,
        extras_require={
            "science": [
                "pymatgen",
                "ase"
                ]
        },
        license='MIT',
        author='Chris Sewell',
        author_email='chrisj_sewell@hotmail.com',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Software Development :: Libraries :: Python Modules',
            'Topic :: Utilities',
        ],
        keywords='python, parser, quantum-espresso, crystal, lammps, gulp',
        zip_safe=True,
        packages=find_packages(),
        package_data={'': ['*.json', '*.csv'],
                      'test_files': ["*"],
                      'schema': ["*"]},
        scripts=['bin/ejplugin_convert']
    )
