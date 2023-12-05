#!/usr/bin/env python
from setuptools import setup, find_packages
setup(
    name="densities-validation",
    author="Blue Brain Project, EPFL",
    setup_requires=["setuptools-git-versioning"],
    setuptools_git_versioning={'enabled': True},
    description=(
        "Performing a set of assertions on cell densities and its subtypes "
        "according to litterature"
    ),
    license="BBP-internal-confidential",
    python_requires=">=3.9",
    install_requires=[
        "voxcell==3.1.6",
    ],
    packages=find_packages(exclude=("tests",)),
    include_package_data=True,
    entry_points={
        "console_scripts": ["densities-validation=validation:main"]
    },
)
