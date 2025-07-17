# GenAIRR: Modular Ig Sequence Simulation

[GenAIRR](https://github.com/MuteJester/GenAIRR) is an AGPL-3 licensed Ig Sequence Simulation framework in Python.

GenAIRR allows users to quickly create complex, modular simulation pipeline using built-in core components (such as Sequence structure, Alleles, Augmenters and Pipelines) or customized implementations; visualize your simulation pipelines, experiment and modify various hyperparameters to control with high precision over the process affecting sequences. Its goal is to allow covering a wide range of scenarios both typical and atypical to aid in the study, benchmarking and development of Ig sequence oriented algorithms such as sequence alignment.

## Features

- Built-in core sequence structure components
- Flexible mutation models and corruption management through "Step" based logic
- New Steps Can be Implemented and Inserted into any Pipeline
- Built-in detailed Data Config files containing empirical distribution and reference data
- Highly sensitive ambiguity resolution framework to ensure reliable ground truth meta information on each simulated sequence

## Documentation Guide

### Getting Started
- **[Step-by-Step Tutorial](step_by_step_tutorial.md)** - Build your first pipeline from scratch
- **[Quick Start Guide](tutorials/Quick Start Guide.ipynb)** - Interactive Jupyter notebook
- **[Getting Started](getting_started.md)** - Basic concepts and examples

### Understanding GenAIRR
- **[Biological Context](biological_context.md)** - What biological processes are being simulated
- **[GenAIRR Flow](genairr_flow.md)** - How the simulation pipeline works
- **[Best Practices](best_practices.md)** - Guidelines for effective use

### Reference Materials  
- **[Parameter Reference](parameter_reference.md)** - Detailed parameter explanations
- **[API Reference](api_reference.md)** - Quick syntax reference
- **[Troubleshooting](troubleshooting.md)** - Common issues and solutions
- **[Migration Guide](migration_guide.md)** - Upgrading from older versions

### Advanced Topics
- **[Advanced Custom Generation](tutorials/Advanced Custom Generation.ipynb)** - Custom allele selection
- **[Introduction to DataConfig](tutorials/Introduction to the DataConfig Object.ipynb)** - Data configuration details
- **[Custom Data Config](custom_data_config.md)** - Using your own germline data

## Using GenAIRR
### Installation Options
To install our latest stable release, run:

```bash
pip install GenAIRR
```
### Resources

For help getting started with Mesa, check out these resources:

- [Getting Started] - Learn about GenAIRR's core concepts and components
- [Advanced Control] - Learn about manipulation of GenAIRR's core components
- [GenAIRR Examples] - Browse various application and example implementations using GenAIRR
- [GitHub Discussions] - Ask questions, make requests and discuss GenAIRR

### Development and Support

GenAIRR is an open source project and welcomes contributions:

- [GitHub Repository] - Access the source code
- [Issue Tracker] - Report bugs or suggest features
- [Contributors Guide] - Learn how to contribute

#### GenAIRR Publication

The original GenAIRR Briefings in Bioinformatics paper is [available here](https://academic.oup.com/bib/article/25/6/bbae556/7863770).


#### Acknowledgments
Some parts of GenAIRR we inspired and adapted from the [AIRRship](https://github.com/Cowanlab/airrship) Package

[contributors guide]: https://github.com/MuteJester/GenAIRR/blob/main/CONTRIBUTING.md
[github repository]: https://github.com/MuteJester/GenAIRR
[github discussions]: https://github.com/MuteJester/GenAIRRdiscussions
[issue tracker]: https://github.com/MuteJester/GenAIRR/issues
[GenAIRR]: https://github.com/MuteJester/GenAIRR
[mesa overview]: overview
[GenAIRR examples]: https://github.com/MuteJester/GenAIRR/tree/master/docs/tutorials
[Getting started]: getting_started
[Advanced Control]: tutorials/Advanced%20Custom%20Generation.ipynb