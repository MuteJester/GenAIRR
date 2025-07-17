# GenAIRR Test Suite Documentation

## Overview
This directory contains a comprehensive test suite for the GenAIRR library, organized into three main test categories that provide thorough coverage of all library functionality.

## Test Structure

### 🧪 Core Functionality Tests (`test_core_functionality.py`)
**Purpose**: Original validation tests ensuring basic library operations work correctly.

**Coverage**:
- Basic sequence simulation and generation
- Pipeline step validation and ordering
- V/D/J position fixing algorithms
- Trimming and corruption functionality
- Mutation rate calculations
- TCR and heavy chain sequence generation

**Test Count**: 23 tests  
**Focus**: Core stability and backward compatibility

### 🔬 Comprehensive Test Suite (`test_comprehensive_suite.py`)
**Purpose**: Enhanced testing covering broader functionality with detailed validation.

**Coverage**:
- Allele creation, properties, and trimming mechanisms
- Mutation models (S5F, Uniform) with various parameters
- Sequence generation variability and quality metrics
- Complete pipeline testing with all augmentation steps
- Data configuration validation and copying
- Utility function testing (translation, weighted choice, stop codon detection)
- Error handling and edge case scenarios
- Light chain vs Heavy chain distinctions
- Integration testing with real-world scenarios
- Performance and memory testing
- Biological constraint validation

**Test Count**: 30 tests  
**Focus**: Comprehensive functionality coverage and quality assurance

### ⚡ Advanced Features Tests (`test_advanced_features.py`)
**Purpose**: Specialized testing for advanced features, performance, and edge cases.

**Coverage**:
- Allele usage distribution analysis
- Configuration serialization/deserialization
- Container data type handling and security
- Custom augmentation step implementation
- Large-scale sequence generation (performance testing)
- Memory efficiency and concurrent execution
- Mutation distribution patterns and statistical analysis
- Pipeline state consistency and dependency validation
- Biological constraint enforcement
- TCR-specific functionality
- Advanced utility function edge cases
- Stress testing and scalability validation

**Test Count**: 21 tests  
**Focus**: Advanced features, performance, and specialized functionality

## Test Execution

### Running All Tests
```bash
python run_all_tests.py
```

### Running Individual Test Suites
```bash
# Core functionality tests
python -m unittest test_core_functionality.TestSequenceSimulation -v

# Comprehensive test suite
python -m unittest test_comprehensive_suite.EnhancedGenAIRRTests -v

# Advanced features tests
python -m unittest test_advanced_features.SpecializedGenAIRRTests -v
```

### Test Runner Features
- **Comprehensive Reporting**: Detailed summary with execution times and success rates
- **Suite Breakdown**: Individual statistics for each test category
- **Error Classification**: Clear distinction between failures and errors
- **Performance Metrics**: Execution timing for performance monitoring
- **Coverage Recommendations**: Suggestions for improving test coverage

## Test Results Interpretation

### Success Rates
- **95-100%**: Excellent - Library is stable and well-tested
- **85-94%**: Good - Minor issues that should be addressed
- **70-84%**: Fair - Significant issues requiring attention
- **Below 70%**: Poor - Major problems that need immediate fixing

### Test Categories
- ✅ **Passed**: Test completed successfully
- ❌ **Failed**: Test assertion failed (expected vs actual mismatch)
- 💥 **Error**: Test encountered an exception during execution
- ⏭️ **Skipped**: Test was skipped (usually due to missing dependencies)

## Development Guidelines

### Adding New Tests
1. **Core Tests**: Add to `test_core_functionality.py` for basic functionality
2. **Comprehensive Tests**: Add to `test_comprehensive_suite.py` for detailed validation
3. **Advanced Tests**: Add to `test_advanced_features.py` for specialized features

### Test Naming Convention
- Use descriptive names that clearly indicate what is being tested
- Follow the pattern: `test_<functionality>_<specific_aspect>`
- Group related tests in logical sections with clear comments

### Test Documentation
- Include comprehensive docstrings explaining test purpose
- Add inline comments for complex test logic
- Document expected outcomes and edge cases

## Dependencies
- **GenAIRR Library**: Main library being tested
- **unittest**: Python standard testing framework
- **numpy**: Numerical operations and array handling
- **tempfile/os**: File system operations for testing
- **threading**: Concurrent execution testing
- **unittest.mock**: Mocking and patching for isolated testing

## Data Files
The test directory includes reference data files:
- `IGHV.fasta`, `IGHD.fasta`, `IGHJ.fasta`, `IGHC.fasta`: Immunoglobulin gene sequences
- `inference_sample.csv`: Sample data for testing inference functionality

## Troubleshooting

### Common Issues
1. **Import Errors**: Ensure GenAIRR library is properly installed
2. **Test Failures**: Check that all dependencies are available
3. **Performance Issues**: Large test batches may require more memory/time

### Environment Setup
```bash
# Ensure virtual environment is activated
source /path/to/venv/bin/activate  # Linux/Mac
# or
\path\to\venv\Scripts\activate     # Windows

# Install required dependencies
pip install -r requirements.txt
```

## Contributing
When contributing new tests:
1. Follow the existing test structure and naming conventions
2. Ensure tests are deterministic and reliable
3. Add appropriate documentation and comments
4. Update this README if adding new test categories
5. Verify all tests pass before submitting

## Contact
For questions about the test suite or to report issues, please refer to the main GenAIRR repository documentation.
