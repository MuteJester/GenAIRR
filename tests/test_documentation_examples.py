#!/usr/bin/env python3
"""
Test script to validate updated documentation examples work correctly.
This script tests the main examples from the updated documentation.
"""

def test_quick_start_example():
    """Test the main Quick Start example from README"""
    print("Testing Quick Start example...")
    
    from GenAIRR.pipeline import AugmentationPipeline
    from GenAIRR.steps import SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity
    from GenAIRR.mutation import S5F
    from GenAIRR.data import HUMAN_IGH_OGRDB
    from GenAIRR.steps.StepBase import AugmentationStep

    # Configure built-in germline data
    AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

    # Build a minimal pipeline
    pipeline = AugmentationPipeline([
        SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
        FixVPositionAfterTrimmingIndexAmbiguity()
    ])

    # Simulate!
    sim = pipeline.execute()
    result = sim.get_dict()
    
    assert 'sequence' in result
    assert 'productive' in result
    assert 'v_call' in result
    print("✓ Quick Start example works correctly")

def test_full_pipeline_example():
    """Test the full heavy chain pipeline example"""
    print("Testing Full Pipeline example...")
    
    from GenAIRR.pipeline import AugmentationPipeline
    from GenAIRR.steps import (
        SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity,
        FixDPositionAfterTrimmingIndexAmbiguity, FixJPositionAfterTrimmingIndexAmbiguity,
        CorrectForVEndCut, CorrectForDTrims, CorruptSequenceBeginning,
        InsertNs, InsertIndels, ShortDValidation, DistillMutationRate
    )
    from GenAIRR.mutation import S5F
    from GenAIRR.data import HUMAN_IGH_OGRDB
    from GenAIRR.steps.StepBase import AugmentationStep

    AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)

    pipeline = AugmentationPipeline([
        SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
        FixVPositionAfterTrimmingIndexAmbiguity(),
        FixDPositionAfterTrimmingIndexAmbiguity(),
        FixJPositionAfterTrimmingIndexAmbiguity(),
        CorrectForVEndCut(),
        CorrectForDTrims(),
        CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
        InsertNs(0.02, 0.5),
        ShortDValidation(),
        InsertIndels(0.5, 5, 0.5, 0.5),
        DistillMutationRate()
    ])
    
    result = pipeline.execute()
    data = result.get_dict()
    
    assert 'sequence' in data
    assert 'mutation_rate' in data
    print("✓ Full Pipeline example works correctly")

def test_naive_sequence_example():
    """Test the naive sequence example"""
    print("Testing Naive Sequence example...")
    
    from GenAIRR.pipeline import AugmentationPipeline
    from GenAIRR.steps import SimulateSequence
    from GenAIRR.mutation import Uniform
    from GenAIRR.data import HUMAN_IGH_OGRDB
    from GenAIRR.steps.StepBase import AugmentationStep

    AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
    
    naive_step = SimulateSequence(Uniform(0, 0), True)
    pipeline = AugmentationPipeline([naive_step])
    naive_seq = pipeline.execute()
    
    result = naive_seq.get_dict()
    assert 'sequence' in result
    assert result['mutation_rate'] == 0.0  # Should be 0 for naive sequence
    print("✓ Naive Sequence example works correctly")

def test_light_chain_example():
    """Test light chain simulation"""
    print("Testing Light Chain example...")
    
    from GenAIRR.pipeline import AugmentationPipeline
    from GenAIRR.steps import (
        SimulateSequence, FixVPositionAfterTrimmingIndexAmbiguity,
        FixJPositionAfterTrimmingIndexAmbiguity, CorrectForVEndCut,
        CorruptSequenceBeginning, InsertNs, InsertIndels, DistillMutationRate
    )
    from GenAIRR.mutation import S5F
    from GenAIRR.data import HUMAN_IGK_OGRDB
    from GenAIRR.steps.StepBase import AugmentationStep

    AugmentationStep.set_dataconfig(HUMAN_IGK_OGRDB)

    pipeline = AugmentationPipeline([
        SimulateSequence(S5F(min_mutation_rate=0.003, max_mutation_rate=0.25), True),
        FixVPositionAfterTrimmingIndexAmbiguity(),
        FixJPositionAfterTrimmingIndexAmbiguity(),
        CorrectForVEndCut(),
        CorruptSequenceBeginning(0.7, [0.4, 0.4, 0.2], 576, 210, 310, 50),
        InsertNs(0.02, 0.5),
        InsertIndels(0.5, 5, 0.5, 0.5),
        DistillMutationRate()
    ])
    
    result = pipeline.execute()
    data = result.get_dict()
    
    assert 'sequence' in data
    assert 'd_call' not in data or data['d_call'] == []  # Light chains don't have D genes
    print("✓ Light Chain example works correctly")

def test_custom_allele_example():
    """Test custom allele selection"""
    print("Testing Custom Allele example...")
    
    from GenAIRR.pipeline import AugmentationPipeline
    from GenAIRR.steps import SimulateSequence
    from GenAIRR.mutation import S5F
    from GenAIRR.data import HUMAN_IGH_OGRDB
    from GenAIRR.steps.StepBase import AugmentationStep

    AugmentationStep.set_dataconfig(HUMAN_IGH_OGRDB)
    
    # Get available alleles (use first available alleles)
    v_family = list(HUMAN_IGH_OGRDB.v_alleles.keys())[0]
    d_family = list(HUMAN_IGH_OGRDB.d_alleles.keys())[0]
    j_family = list(HUMAN_IGH_OGRDB.j_alleles.keys())[0]
    
    custom_step = SimulateSequence(
        S5F(0.003, 0.25),
        True,
        specific_v=HUMAN_IGH_OGRDB.v_alleles[v_family][0],
        specific_d=HUMAN_IGH_OGRDB.d_alleles[d_family][0],
        specific_j=HUMAN_IGH_OGRDB.j_alleles[j_family][0]
    )
    
    pipeline = AugmentationPipeline([custom_step])
    result = pipeline.execute()
    data = result.get_dict()
    
    assert 'sequence' in data
    assert 'v_call' in data
    print("✓ Custom Allele example works correctly")

if __name__ == "__main__":
    print("Testing updated documentation examples...")
    print("=" * 50)
    
    try:
        test_quick_start_example()
        test_full_pipeline_example()
        test_naive_sequence_example()
        test_light_chain_example()
        test_custom_allele_example()
        
        print("=" * 50)
        print("✅ All documentation examples work correctly!")
        print("The updated documentation syntax is validated.")
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
