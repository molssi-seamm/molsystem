#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy>=1.16.4',
    'pandas>=0.24.2',
]

setup_requirements = [
    'pytest-runner>=5.1',
]

test_requirements = [
    'pytest>=4.5.0',
]

setup(
    name='molsystem',
    version='0.1.1',
    description=("Molsystem provides a general class for handling "
                 "molecular and periodic systems"),
    long_description=readme + '\n\n' + history,
    author="Paul Saxe",
    author_email='psaxe@vt.edu',
    url='https://github.com/paulsaxe/molsystem',
    packages=find_packages(include=['molsystem']),
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='molsystem',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
