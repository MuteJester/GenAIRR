from enum import Enum, auto
from dataclasses import dataclass
from GenAIRR.data import builtin_heavy_chain_data_config, builtin_kappa_chain_data_config, builtin_lambda_chain_data_config, builtin_tcrb_data_config
from GenAIRR.dataconfig import DataConfig


class ChainType(Enum):
    BCR_HEAVY = auto()
    BCR_LIGHT_KAPPA = auto()
    BCR_LIGHT_LAMBDA = auto()
    #TCR_ALPHA = auto()
    TCR_BETA = auto()
    # TCR_GAMMA = auto()
    # TCR_DELTA = auto()

@dataclass(frozen=True)
class ChainInfo:
    name: str
    has_d: bool
    dataconfig: DataConfig

# Metadata registry for each chain type
CHAIN_TYPE_INFO = {
    ChainType.BCR_HEAVY: ChainInfo(name="BCR Heavy Chain", has_d=True,dataconfig=builtin_heavy_chain_data_config()),
    ChainType.BCR_LIGHT_KAPPA: ChainInfo(name="BCR Light Chain Kappa", has_d=False,dataconfig=builtin_kappa_chain_data_config()),
    ChainType.BCR_LIGHT_LAMBDA: ChainInfo(name="BCR Light Chain Lambda", has_d=False,dataconfig=builtin_lambda_chain_data_config()),
    #ChainType.TCR_ALPHA: ChainInfo(name="TCR Alpha Chain", has_d=False),
    ChainType.TCR_BETA: ChainInfo(name="TCR Beta Chain", has_d=True,dataconfig=builtin_tcrb_data_config()),
    # ChainType.TCR_GAMMA: ChainInfo(name="TCR Gamma Chain", has_d=False),
    # ChainType.TCR_DELTA: ChainInfo(name="TCR Delta Chain", has_d=True),
}
