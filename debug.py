import pickle
import os
import sys

# STEP 1: Patch the module path (fake the old one)
import GenAIRR.dataconfig as new_module

sys.modules["GenAIRR.utilities.data_config"] = new_module

# STEP 2: Paths to files
data_dir = r"C:\Users\tomas\Desktop\AlignAIRR\tests"
file_map = {
    'Genotyped_DataConfig.pkl': 'Genotyped_DataConfig.pkl',
}

# STEP 3: Load and re-save
for filename in file_map:
    full_path = os.path.join(data_dir, filename)

    # Load with old module path alias
    with open(full_path, 'rb') as f:
        obj = pickle.load(f)

    # Re-save with updated module path
    with open(full_path, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Resaved: {filename}")
