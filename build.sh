#!/bin/bash

# Create a new gpatch release and update all repositories. Takes one
# argument to indicate whether this is a major, minor, or patch release.

# Increment version numbers
bump2version --list $1

# Update PyPI
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/patch_genome-$(cat VERSION)*

# Update BioConda
