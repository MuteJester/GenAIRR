from setuptools import setup, find_packages

setup(
    name='GenAIRR',
    version='0.1.0',
    author='Thomas Konstantinovsky & Ayelet Peres',
    author_email='thomaskon90@gmail.com',
    description='Advanced IG Sequence Simulation Suit for Alignment Model Benchmarking',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/MuteJester/GenAIRR',
    package_dir={'': 'src'},  # Tells setuptools to look for packages in the src directory
    packages=find_packages(where='src'),  # Automatically find and include all packages under src
    install_requires=[
        # List your project's dependencies here
        # e.g., 'numpy>=1.18.0', 'pandas>=1.0.0'
    ],
    classifiers=[
        # Choose your license as you wish
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.7',
)
