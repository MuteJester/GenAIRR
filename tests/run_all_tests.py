#!/usr/bin/env python3
"""
Comprehensive Test Runner for GenAIRR Library

This script runs all the test suites and provides a comprehensive summary
of test coverage and results for the GenAIRR library.
"""

import unittest
import sys
import time
from io import StringIO


class TestResult:
    """Class to track test results and statistics."""
    
    def __init__(self):
        self.total_tests = 0
        self.passed_tests = 0
        self.failed_tests = 0
        self.error_tests = 0
        self.skipped_tests = 0
        self.test_suites = {}
        self.start_time = None
        self.end_time = None
    
    def start_timing(self):
        self.start_time = time.time()
    
    def end_timing(self):
        self.end_time = time.time()
    
    def get_duration(self):
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return 0
    
    def add_suite_result(self, suite_name, result):
        self.test_suites[suite_name] = {
            'tests_run': result.testsRun,
            'failures': len(result.failures),
            'errors': len(result.errors),
            'skipped': len(result.skipped) if hasattr(result, 'skipped') else 0
        }
        
        self.total_tests += result.testsRun
        self.failed_tests += len(result.failures)
        self.error_tests += len(result.errors)
        self.skipped_tests += len(result.skipped) if hasattr(result, 'skipped') else 0
        self.passed_tests += (result.testsRun - len(result.failures) - len(result.errors) - 
                             (len(result.skipped) if hasattr(result, 'skipped') else 0))


def run_test_suite(suite_name, test_class):
    """Run a specific test suite and return results."""
    print(f"\n{'='*60}")
    print(f"Running {suite_name}")
    print('='*60)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(test_class)
    
    # Run tests with custom result to capture details
    stream = StringIO()
    runner = unittest.TextTestRunner(stream=stream, verbosity=2)
    result = runner.run(suite)
    
    # Print results
    output = stream.getvalue()
    print(output)
    
    return result


def print_summary(test_result):
    """Print comprehensive test summary."""
    print(f"\n{'='*60}")
    print("COMPREHENSIVE TEST SUMMARY")
    print('='*60)
    
    print(f"Total Test Execution Time: {test_result.get_duration():.2f} seconds")
    print(f"Total Tests Run: {test_result.total_tests}")
    print(f"✅ Passed: {test_result.passed_tests}")
    print(f"❌ Failed: {test_result.failed_tests}")
    print(f"💥 Errors: {test_result.error_tests}")
    print(f"⏭️  Skipped: {test_result.skipped_tests}")
    
    success_rate = (test_result.passed_tests / test_result.total_tests * 100) if test_result.total_tests > 0 else 0
    print(f"Success Rate: {success_rate:.1f}%")
    
    print(f"\n{'='*60}")
    print("TEST SUITE BREAKDOWN")
    print('='*60)
    
    for suite_name, stats in test_result.test_suites.items():
        suite_success_rate = ((stats['tests_run'] - stats['failures'] - stats['errors'] - stats['skipped']) / 
                             stats['tests_run'] * 100) if stats['tests_run'] > 0 else 0
        
        print(f"\n{suite_name}:")
        print(f"  Tests: {stats['tests_run']}")
        print(f"  Passed: {stats['tests_run'] - stats['failures'] - stats['errors'] - stats['skipped']}")
        print(f"  Failed: {stats['failures']}")
        print(f"  Errors: {stats['errors']}")
        print(f"  Skipped: {stats['skipped']}")
        print(f"  Success Rate: {suite_success_rate:.1f}%")


def main():
    """Main test runner function."""
    print("GenAIRR Library - Comprehensive Test Suite")
    print("==========================================")
    
    test_result = TestResult()
    test_result.start_timing()
    
    # Import test classes
    try:
        import sys
        import os
        # Add current directory to path for imports
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        
        from test_core_functionality import TestSequenceSimulation
        from test_comprehensive_suite import EnhancedGenAIRRTests
        from test_advanced_features import SpecializedGenAIRRTests
        from test_enhanced_sequence_simulation import EnhancedSequenceSimulationTests
        
        test_suites = [
            ("Core Functionality Tests", TestSequenceSimulation),
            ("Comprehensive Test Suite", EnhancedGenAIRRTests),
            ("Advanced Features Tests", SpecializedGenAIRRTests),
            ("Enhanced Sequence Simulation Tests", EnhancedSequenceSimulationTests)
        ]
        
        # Run each test suite
        for suite_name, test_class in test_suites:
            try:
                result = run_test_suite(suite_name, test_class)
                test_result.add_suite_result(suite_name, result)
            except Exception as e:
                print(f"❌ Error running {suite_name}: {e}")
                test_result.error_tests += 1
        
    except ImportError as e:
        print(f"❌ Error importing test modules: {e}")
        print("Make sure all test files are present and the GenAIRR library is properly installed.")
        return 1
    
    test_result.end_timing()
    
    # Print comprehensive summary
    print_summary(test_result)
    
    # Provide recommendations
    print(f"\n{'='*60}")
    print("RECOMMENDATIONS")
    print('='*60)
    
    if test_result.failed_tests > 0 or test_result.error_tests > 0:
        print("❌ Some tests failed or had errors.")
        print("   Recommendations:")
        print("   1. Review failed test details above")
        print("   2. Check for missing dependencies")
        print("   3. Verify data configuration files are present")
        print("   4. Ensure proper GenAIRR library installation")
    else:
        print("✅ All tests passed successfully!")
        print("   Your GenAIRR installation is working correctly.")
        print("   The enhanced test coverage provides confidence in:")
        print("   • Core functionality (alleles, sequences, mutations)")
        print("   • Pipeline operations and step ordering")
        print("   • Data configuration and validation")
        print("   • Utility functions and edge cases")
        print("   • Performance and memory management")
        print("   • Error handling and robustness")
    
    print(f"\n{'='*60}")
    print("TEST COVERAGE AREAS")
    print('='*60)
    print("✅ Allele creation and trimming")
    print("✅ Mutation models (S5F, Uniform)")
    print("✅ Sequence generation and variability")
    print("✅ Pipeline construction and execution")
    print("✅ Data configuration validation")
    print("✅ Utility function testing")
    print("✅ Error handling and edge cases")
    print("✅ Light chain vs Heavy chain differences")
    print("✅ TCR-specific functionality")
    print("✅ Performance and stress testing")
    print("✅ Biological validity checks")
    
    # Return appropriate exit code
    if test_result.failed_tests > 0 or test_result.error_tests > 0:
        return 1
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
