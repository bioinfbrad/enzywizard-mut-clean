#!/usr/bin/env python
from setuptools import setup, find_packages
import os

# Read the version from version.py without importing the package
version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_mut_clean', 'version.py')
with open(version_file) as f:
    exec(f.read())  # defines __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-mut-clean",
    version=__version__,                     # dynamically read from version.py (1.0.1)
    author="bioinfbrad",
    description=(
        "Clean a wild‑type protein structure and its corresponding mutant structure, "
        "remap amino acid substitutions, generate multi‑format cleaned files (CIF, PDB, FASTA), "
        "and produce a detailed traceable JSON report."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-mut-clean",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.86",          # For structure I/O and sequence handling
        "numpy>=1.23.5",            # Numerical backend (used in clean_algorithms)
        "openmm>=8.5.0",            # Molecular mechanics engine (hydrogen addition, force fields)
        "pdbfixer>=1.12",           # PDBFixer APIs for cleaning, heterogen removal, missing atom repair
        # 'packaging' and 'scipy' are not directly used by this package,
        # but may be pulled in as transitive dependencies via Conda.
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-mut-clean = enzywizard_mut_clean.cli:main",
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
