from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='GenAIRR',
    version='0.6.0',
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
            'data/mutation_model_parameters/*.json'
        ]
    },  # Include data files from subfolders
    include_package_data=True,  # Include everything in source control
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='immunogenetics, sequence simulation, bioinformatics, alignment benchmarking',
    python_requires='>=3.9,<3.13',
)
