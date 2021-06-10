import setuptools

#! /usr/bin/python
# -*- coding: utf-8 -*-

from typing import List
from setuptools import setup, find_packages

PACKAGE_VERSION = '1.3.2'

def write_version_py(filename='polygenic/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SETUP.PY
__version__ = '{}'
"""
    with open(filename, 'w') as f:
        f.write(cnt.format(PACKAGE_VERSION))


write_version_py()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="polygenic",
    version=PACKAGE_VERSION,
    author="Marcin Piechota, Wojciech Gałan",
    author_email="piechota@intelliseq.com",
    description="Polygenic score computation",
    #long_description=long_description,
    #long_description_content_type="text/reStructuredText",
    url="https://github.com/marpiech/polygenic",
    packages=setuptools.find_packages(),
    package_data={'polygenic': ['*.cfg']},
    license="Intelliseq dual licenses this package. For commercial use, please contact [contact @ intelliseq.com](mailto:contact@intelliseq.com). For non-commercial use, this license permits use of the software only by government agencies, schools, universities, non-profit organizations or individuals on projects that do not receive external funding other than government research grants and contracts. Any other use requires a commercial license. For the full license, please see [LICENSE.md](https://github.com/intelliseq/polygenic/blob/master/LICENSE.md), in this source repository.",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: Free for non-commercial use',
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    install_requires=[
        'namedlist==1.8',
        'numpy==1.20.3'
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'polygenic=polygenic.polygenic:main',
        ],
    },
    test_suite='nose.collector',
    tests_require=['nose>=1.0'],
    setup_requires=['nose>=1.0'],
)
