import os
import subprocess
from pathlib import Path

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


# ── CMake extension build ────────────────────────────────────────

class CMakeExtension(Extension):
    """A dummy extension that triggers the CMake build."""

    def __init__(self, name, cmake_source_dir=""):
        super().__init__(name, sources=[])
        self.cmake_source_dir = os.path.abspath(cmake_source_dir)


class CMakeBuild(build_ext):
    """Custom build_ext that drives CMake to build the C shared library."""

    def run(self):
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake is required to build GenAIRR. "
                "Install it with: pip install cmake"
            )
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        # Determine where the built library should go
        ext_fullpath = Path(self.get_ext_fullpath(ext.name)).resolve()
        output_dir = ext_fullpath.parent

        build_dir = Path(self.build_temp).resolve() / "cmake_build"
        build_dir.mkdir(parents=True, exist_ok=True)

        cfg = "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={output_dir}",
            f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={output_dir}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DGENAIRR_BUILD_TESTS=OFF",
            "-DGENAIRR_BUILD_BENCHMARKS=OFF",
        ]

        build_args = ["--config", cfg]

        # Parallel build
        if hasattr(os, "cpu_count"):
            n_jobs = os.cpu_count() or 1
            build_args += [f"-j{n_jobs}"]

        subprocess.check_call(
            ["cmake", ext.cmake_source_dir] + cmake_args,
            cwd=build_dir,
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args,
            cwd=build_dir,
        )


# ── Package metadata ────────────────────────────────────────────

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='GenAIRR',
    version='1.0.0',
    author='Thomas Konstantinovsky & Ayelet Peres',
    author_email='thomaskon90@gmail.com',
    description='An advanced immunoglobulin sequence simulation suite for benchmarking alignment models and sequence analysis.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/MuteJester/GenAIRR',
    project_urls={
        "Bug Tracker": "https://github.com/MuteJester/GenAIRR/issues"
    },
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    package_data={
        'GenAIRR': [
            'data/*.pkl',
            'data/*.json',
            'data/builtin_dataconfigs/*.pkl',
            'data/builtin_dataconfigs/*.json',
            'data/mutation_model_parameters/*.pkl',
            'data/mutation_model_parameters/*.json',
        ]
    },
    include_package_data=True,
    install_requires=[],
    ext_modules=[
        CMakeExtension(
            "GenAIRR._native._cmake",
            cmake_source_dir="src/GenAIRR/_native/csrc",
        ),
    ],
    cmdclass={"build_ext": CMakeBuild},
    entry_points={
        'console_scripts': [
            'genairr-mcp=GenAIRR.mcp_server:main',
        ],
    },
    extras_require={
        'dataconfig': ['numpy>=1.24', 'scipy>=1.10'],
        'viz': ['graphviz>=0.20'],
        'mcp': ['fastmcp>=2.0'],
        'all': ['numpy>=1.24', 'scipy>=1.10', 'graphviz>=0.20', 'tqdm>=4.60', 'fastmcp>=2.0'],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: C',
    ],
    keywords='immunogenetics, sequence simulation, bioinformatics, alignment benchmarking',
    python_requires='>=3.9',
)
