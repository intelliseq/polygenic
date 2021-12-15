#!/bin/bash
python3 setup.py sdist bdist_wheel
VERSION=$(cat polygenic/version.py | grep version | grep -oP "[0-9]+\.[0-9]+\.[0-9]+")
twine upload --repository pypi dist/polygenic-$VERSION*
