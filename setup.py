"""GenAIRR — pure-Python package metadata.

The simulation kernel lives in the separate ``genairr_engine`` Rust
crate at ``engine_rs/``, built and distributed by maturin. This
``setup.py`` only packages the Python wrappers + reference data.
"""
from setuptools import find_packages, setup


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="GenAIRR",
    version="1.0.0",
    author="Thomas Konstantinovsky & Ayelet Peres",
    author_email="thomaskon90@gmail.com",
    description=(
        "Synthetic immune-receptor-sequence simulator for benchmarking "
        "alignment models and sequence analysis."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MuteJester/GenAIRR",
    project_urls={
        "Bug Tracker": "https://github.com/MuteJester/GenAIRR/issues",
    },
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={
        "GenAIRR": [
            "data/*.pkl",
            "data/*.json",
            "data/builtin_dataconfigs/*.pkl",
            "data/builtin_dataconfigs/*.json",
            "data/mutation_model_parameters/*.pkl",
            "data/mutation_model_parameters/*.json",
        ],
    },
    include_package_data=True,
    # ``genairr_engine`` is a separate maturin-built wheel. Pin it
    # at install time so users don't have to remember.
    install_requires=[
        "genairr_engine>=0.0.1",
    ],
    extras_require={
        "dataconfig": ["numpy>=1.24", "scipy>=1.10"],
        "viz": ["graphviz>=0.20"],
        "all": [
            "numpy>=1.24",
            "scipy>=1.10",
            "graphviz>=0.20",
            "tqdm>=4.60",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Rust",
    ],
    keywords="immunogenetics, sequence simulation, bioinformatics, alignment benchmarking",
    python_requires=">=3.9",
)
