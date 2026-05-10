---
title: Help & Troubleshooting
sidebar_label: Overview
---

# Help & Troubleshooting

This page provides answers to common questions and solutions for issues you might encounter while using GenAIRR.

## Frequently Asked Questions

### Why is `mutate()` failing on my TCR experiment?
Somatic hypermutation (SHM) is a biological phenomenon specific to B-cells. GenAIRR enforces this biological constraint by raising a `ValueError` if you attempt to call `.mutate()` on an experiment configured for T-cell receptors (like `human_trb`).

### How can I generate only productive sequences?
GenAIRR uses **Contracts** to ensure sequences meet productivity requirements by construction. Use the `productive()` contract bundle when running your experiment:
```python
import GenAIRR as ga
result = exp.run(n=100, respect=ga.productive())
```

### Are simulations reproducible across different computers?
Yes. GenAIRR is bit-for-bit deterministic. Using the same **Seed**, **Reference Data**, and **Pipeline** will produce identical results on Linux, macOS, and Windows.

## Common Issues

### `ModuleNotFoundError: No module named 'GenAIRR._engine'`
This error indicates that the core simulation kernel is not correctly installed or accessible in your current Python environment. This can happen if the installation was interrupted or if you are using an incompatible Python version.

**Solution:** Try forcing a reinstallation of the package:
```bash
pip install GenAIRR --force-reinstall
```

### Simulation performance is lower than expected
GenAIRR is designed to produce tens of thousands of sequences per second. If you are experiencing significantly lower performance, ensure that:
1. You are not running in a debugger or with heavy profiling enabled.
2. Your system has sufficient memory available for the batch size you are simulating.
3. You are using the official pre-built wheels (which are optimized for your platform).

## Getting Further Help
If you encounter a bug or have a feature request, please open an issue on our [GitHub repository](https://github.com/im-repertoire/GenAIRR/issues).
