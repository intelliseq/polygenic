from typing import Dict
from typing import List
from typing import Union
from typing import NamedTuple


class SnpData(NamedTuple):
    ref:str
    alts:List[str]
    genotype:List[str]


class SnpDataManySamples(NamedTuple):
    ref: str
    alts: List[str]
    genotypes: Dict[str, List[str]] # {sample: [A, G]}


# class ModelDescriptionInfo(NamedTuple):
#     model_fname:str
#     model_path: str
#     desc_paths:List[str]
#     plot_data_path:Union[str, None] = None