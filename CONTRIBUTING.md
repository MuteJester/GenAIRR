# Contributing to GenAIRR

As an open-source project, **GenAIRR** welcomes contributions in various forms, catering to all levels of expertise, from beginners to experienced developers. Our goal is to foster a collaborative environment where everyone can contribute and learn. Whether you’re interested in improving the documentation, submitting bug reports, or adding new features, your input is valuable to us.

## Contribution Opportunities

Contributions can come in many forms, including:

- Code patches
- Bug reports and reviews of patches
- New feature development
- Enhancements to documentation
- Tutorials and examples

No contribution is too small! However, large changes should be discussed beforehand via our [GitHub Discussions] or by opening [an issue].

## Release Process

GenAIRR uses automated PyPI publishing through GitHub Actions. When you create a new release, the package will be automatically published to PyPI.

### Creating a Release

You can trigger automatic PyPI publishing in two ways:

#### Option 1: Create a GitHub Release (Recommended)
1. Go to the [Releases page](https://github.com/MuteJester/GenAIRR/releases)
2. Click "Create a new release"
3. Create a new tag with the version number (e.g., `v0.6.0`)
4. Add release notes describing the changes
5. Click "Publish release"

#### Option 2: Push a Git Tag
1. Update the version in `setup.py`
2. Commit your changes
3. Create and push a git tag:
   ```bash
   git tag v0.6.0
   git push origin v0.6.0
   ```

### Prerequisites

For automatic PyPI publishing to work, the repository must have:
- A `PYPI_API_TOKEN` secret configured in the GitHub repository settings
- The token should have permission to publish to the GenAIRR package on PyPI

### What Happens Automatically

When a release is created or a tag is pushed:
1. GitHub Actions will trigger the `Upload Python Package` workflow
2. The workflow will:
   - Set up a Python environment
   - Install build dependencies
   - Build the package (both source distribution and wheel)
   - Publish to PyPI using the stored API token

### Troubleshooting

If the automatic publishing fails:
1. Check the Actions tab for error details
2. Verify the `PYPI_API_TOKEN` secret is correctly configured
3. Ensure the version number in `setup.py` matches the tag/release
4. Check that the version doesn't already exist on PyPI
