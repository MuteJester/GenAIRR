#!/usr/bin/env python3
"""
GenAIRR Test Configuration

This module contains configuration settings and utilities for the GenAIRR test suite.
"""

# Test suite configuration
TEST_SUITES = {
    'core': {
        'module': 'test_core_functionality',
        'class': 'TestSequenceSimulation',
        'description': 'Core Functionality Tests',
        'focus': 'Basic operations and backward compatibility'
    },
    'comprehensive': {
        'module': 'test_comprehensive_suite',
        'class': 'EnhancedGenAIRRTests',
        'description': 'Comprehensive Test Suite',
        'focus': 'Detailed functionality coverage and quality assurance'
    },
    'advanced': {
        'module': 'test_advanced_features',
        'class': 'SpecializedGenAIRRTests',
        'description': 'Advanced Features Tests',
        'focus': 'Performance, edge cases, and specialized functionality'
    }
}

# Test execution settings
DEFAULT_VERBOSITY = 2
PERFORMANCE_TEST_TIMEOUT = 300  # seconds
LARGE_BATCH_SIZE = 1000
STRESS_TEST_ITERATIONS = 500

# Test data configuration
TEST_DATA_FILES = [
    'IGHV.fasta',
    'IGHD.fasta', 
    'IGHJ.fasta',
    'IGHC.fasta',
    'inference_sample.csv'
]

# Expected success rates
EXPECTED_SUCCESS_RATES = {
    'core': 95.0,
    'comprehensive': 100.0,
    'advanced': 95.0,
    'overall': 95.0
}
