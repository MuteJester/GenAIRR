import pickle
import os

module_dir = os.path.dirname(__file__)

def builtin_heavy_chain_data_config():
    data_path = os.path.join(module_dir, 'HeavyChain_DataConfig_OGRDB_V3.pkl')
    with open(data_path, 'rb') as h:
        return pickle.load(h)

def builtin_kappa_chain_data_config():
    data_path = os.path.join(module_dir, 'LightChain_KAPPA_DataConfigV3.pkl')
    with open(data_path, 'rb') as h:
        return pickle.load(h)

def builtin_lambda_chain_data_config():
    data_path = os.path.join(module_dir, 'LightChain_LAMBDA_DataConfigV3.pkl')
    with open(data_path, 'rb') as h:
        return pickle.load(h)