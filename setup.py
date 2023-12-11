import os
from setuptools import setup, find_packages

COMMIT_TAG = os.environ.get("CI_COMMIT_TAG")

setup(
    name="densities-validation",
    author="Blue Brain Project, EPFL",
    setup_requires=["setuptools_scm"],
    use_scm_version=True if not COMMIT_TAG else None,
    version=COMMIT_TAG if COMMIT_TAG else None,
    description=(
        "Performing a set of assertions on cell densities and its subtypes "
        "according to litterature"
    ),
    license="BBP-internal-confidential",
    python_requires=">=3.8",
    install_requires=[
        "voxcell==3.1.6",
    ],
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": ["densities-validation=src.validation:main"]
    },
)
