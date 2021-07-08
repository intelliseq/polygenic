import math
import logging
import operator
import functools
import numpy as np
from typing import Set
from typing import Callable
from typing import Iterable
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from typing import NamedTuple
from polygenic.lib.data_access.allele_frequency_accessor import AlleleFrequencyAccessor
from polygenic.seqql.category import validate_categories_names
from polygenic.seqql.category import validate_categories_ranges
from polygenic.seqql.category import QuantitativeCategory
from polygenic.lib.data_access.data_accessor import VcfAccessor
from polygenic.lib.data_access.dto import SnpData
from polygenic.lib.data_access.data_accessor import DataNotPresentError

logger = logging.getLogger('description_language.' + __name__)


class FunctionData(NamedTuple):
    function_to_apply_on_allele_pair: Callable
    function_to_apply_between_allele_pairs: Callable
    neutral_number: float
    function_to_compute_ranges: Callable


def beta(alleles: Iterable[str], effect_allele: str, beta_value: float) -> float:
    return alleles.count(effect_allele) * beta_value


def odds_ratio(alleles: Iterable[str], effect_allele: str, odds_ratio_value: float) -> float:
    if alleles[0] == alleles[1] == effect_allele:
        return odds_ratio_value
    elif effect_allele in alleles:
        return 1
    else:
        return 1 / odds_ratio_value


def prbs(alleles: Iterable[str], effect_allele: str, odds_ratio_value: float) -> float:
    return alleles.count(effect_allele) * math.log10(odds_ratio_value)


def sub(a, b):
    return a - b


def div_and_ln(a, b):
    return math.log(a / b)


FUNCTION_DICT = {
    'beta': FunctionData(beta, sum, 0, sub),
    'odds_ratio': FunctionData(odds_ratio, functools.partial(functools.reduce, operator.mul), 1, div_and_ln),
    'prbs': FunctionData(prbs, sum, 0, sub)
}


class ConflictingAlleleBetweenDataAndModel(RuntimeError):
    pass


class ModelData(NamedTuple):
    effect_allele: str
    coeff_value: float


class PolygenicRiskScoreResult(NamedTuple):
    initial_range_fraction: float
    imputed_range_fraction: float
    pop_freq_range_fraction: float
    missing_data_range_fraction: float
    initial_result: float
    imputed_result: float
    pop_freq_result: float
    missing_data_result: float
    result: float
    scaled_result: float
    category: str
    boundaries: Dict[str, Dict[str, float]]


class Data(object):
    def __init__(self, vcf_data_accessor: VcfAccessor,
                 allele_freq_in_population_accessor: VcfAccessor, sample_name: str, population: str, model):
        self.__vcf_data_accessor = vcf_data_accessor
        self.__allele_freq_in_population_accessor = allele_freq_in_population_accessor
        self._sample_name = sample_name
        self._population = population
        self._initial_data: Dict[str, SnpData] = {}
        self._imputed_data: Dict[str, SnpData] = {}
        self._allele_freq_in_population: Dict[str, Tuple[str]] = {}
        self._rsids_absent_in_allele_freq_data: Set[str] = set()
        self._model = model
        self._read_data()

    @property
    def sample_name(self):
        return self._sample_name

    @property
    def population(self):
        return self._population

    @property
    def initial_data(self):
        return self._initial_data

    @property
    def imputed_data(self):
        return self._imputed_data

    @property
    def allele_freq_in_population(self):
        return self._allele_freq_in_population

    @property
    def rsids_absent_in_allele_freq_data(self):
        return self._rsids_absent_in_allele_freq_data

    def _read_data(self) -> None:
        required_rsids = self._model.get_rsids()
        #print(required_rsids)
        logger.debug('Looking for rsids in initial data')
        self._initial_data = rsids_data_from_data_source(required_rsids, self.__vcf_data_accessor,
                                                         self._sample_name, False)
        self._model.validate_allele_collection(self._initial_data, 'initial data')
        rsids_not_found_in_initial_data: Set[str] = required_rsids - set(self._initial_data)
        #print(rsids_not_found_in_initial_data)
        logger.debug('Looking for rsids in initial data')
        self._imputed_data = rsids_data_from_data_source(rsids_not_found_in_initial_data, self.__vcf_data_accessor,
                                                         self._sample_name, True)
        self._model.validate_allele_collection(self._imputed_data, 'imputed data')
        rsids_not_found_in_imputed_data: Set[str] = rsids_not_found_in_initial_data - set(self._imputed_data)
        #print(rsids_not_found_in_imputed_data)
        logger.debug('Looking for rsids in allele frequency data')
        for rsid in rsids_not_found_in_imputed_data:
            #print(rsid)
            try:
                #print("====> AF " + rsid + " " + self._population)
                #print(self.__allele_freq_in_population_accessor)
                allele_frequencies: Dict[str, float] = self.__allele_freq_in_population_accessor.get_af_by_pop(rsid, self._population)
                #print(allele_frequencies)
                alleles = set(allele_frequencies)
                validate_alleles(self._model.snips_and_coefficients[rsid].effect_allele, alleles, rsid,
                                 'allele frequency data')
                self._allele_freq_in_population[rsid] = max(hardy_weinberg_diploid_frequencies(allele_frequencies).items(), key=lambda x: x[1])[0]
                
            except KeyError:
                #print("keyerror")
                logger.debug(
                    f'Either rsid: {rsid} or population: {self._population} not found in allele frequency data')
                self._rsids_absent_in_allele_freq_data.add(rsid)
                continue
            except AttributeError:
                #print("attributeerror")
                logger.debug(
                    f'Either rsid: {rsid} or population: {self._population} not found in allele frequency data')
                self._rsids_absent_in_allele_freq_data.add(rsid)
                continue
            except DataNotPresentError:
                #print("datanotpresenterror")
                logger.debug(
                    f'Either rsid: {rsid} or population: {self._population} not found in allele frequency data')
                self._rsids_absent_in_allele_freq_data.add(rsid)
                continue


    def compute_model(self) -> PolygenicRiskScoreResult:

        function_data = FUNCTION_DICT[self._model.model_type]
        result_for_initial_data = compute_part_of_the_model(self._model.snips_and_coefficients, self.initial_data,
                                                            function_data)
        range_for_initial = compute_range(
            (self._model.snips_and_coefficients[rsid] for rsid, _ in self.initial_data.items()), function_data)
        #print(self.initial_data.items().__len__())
        #print(self.imputed_data.items().__len__())
        #print(self.allele_freq_in_population.items().__len__())
        #print(self.rsids_absent_in_allele_freq_data.__len__())
        result_for_imputed_data = compute_part_of_the_model(self._model.snips_and_coefficients, self.imputed_data,
                                                            function_data)
        range_for_imputed = compute_range((
            self._model.snips_and_coefficients[rsid] for rsid, _ in self.imputed_data.items()), function_data)
        result_for_allele_freq_data = compute_part_of_the_model(self._model.snips_and_coefficients,
                                                                self.allele_freq_in_population,
                                                                function_data)
        range_for_allele_freq = compute_range((
            self._model.snips_and_coefficients[rsid] for rsid, _ in self.allele_freq_in_population.items()),
            function_data)
        result_for_not_found = function_data.neutral_number
        range_for_not_found = compute_range((
            self._model.snips_and_coefficients[rsid] for rsid in self.rsids_absent_in_allele_freq_data), function_data)

        range_for_model = sum([range_for_initial, range_for_imputed, range_for_allele_freq, range_for_not_found])
        result_for_model = function_data.function_to_apply_between_allele_pairs(
            [result_for_initial_data, result_for_imputed_data, result_for_allele_freq_data, result_for_not_found])
        category_boundaries = self._model.get_category_boundaries()
        # category_boundaries = transform_boundaries_if_needed(self._model, self._model.get_category_boundaries())

        return PolygenicRiskScoreResult(
            range_for_initial / range_for_model,
            range_for_imputed / range_for_model,
            range_for_allele_freq / range_for_model,
            range_for_not_found / range_for_model,
            result_for_initial_data,
            result_for_imputed_data,
            result_for_allele_freq_data,
            result_for_not_found,
            result_for_model,
            get_category(self._model.categories, result_for_model).scale_cat(result_for_model),
            get_category(self._model.categories, result_for_model).name,
            category_boundaries
        )


class PolygenicRiskScore(object):
    def __init__(self, categories: Iterable[QuantitativeCategory], snips_and_coefficients: Dict[str, ModelData],
                 model_type: str):
        super().__init__()
        implemented_model_types = {'beta', 'odds_ratio', 'prbs'}
        self.model_type = validate_model_type(model_type, implemented_model_types)
        self.categories = categories
        self.__validate_categories()
        self.snips_and_coefficients = snips_and_coefficients
        self.__data_already_read = False

    def get_rsids(self) -> Set[str]:
        return set(self.snips_and_coefficients.keys())

    def __validate_categories(self):
        validate_categories_names(self.categories)
        validate_categories_ranges(self.categories)

    def validate_allele_collection(self, allele_collection: Dict[str, SnpData], data_source: str):
        for rsid, snp_data in self.snips_and_coefficients.items():
            if rsid in allele_collection:
                alleles_in_data = set(allele_collection[rsid].alts)
                alleles_in_data.add(allele_collection[rsid].ref)
                validate_alleles(snp_data.effect_allele, alleles_in_data, rsid, data_source)

    def get_category_boundaries(self) -> Dict[str, Dict[str, float]]:
        return {category.name: {'from':category.from_, 'to': category.to} for category in self.categories}


def get_category(categories: List[QuantitativeCategory], num: float):
    for category in categories:
        if category.from_ <= num <= category.to:
            return category
    raise RuntimeError("No category matched")


def compute_boundaries(model_datas: Iterable[ModelData], function_data: FunctionData) -> Tuple[float]:
    positive_contributions = []
    negative_contributions = []
    for model_data in model_datas:
        extreme_results = [
            function_data.function_to_apply_on_allele_pair(2 * ['.'], model_data.effect_allele, model_data.coeff_value),
            function_data.function_to_apply_on_allele_pair(2 * [model_data.effect_allele], model_data.effect_allele,
                                                           model_data.coeff_value)
        ]
        extreme_results.sort()
        negative_contributions.append(extreme_results[0])
        positive_contributions.append(extreme_results[1])
    if negative_contributions:
        return function_data.function_to_apply_between_allele_pairs(
            negative_contributions), function_data.function_to_apply_between_allele_pairs(positive_contributions)
    else:
        return function_data.neutral_number, function_data.neutral_number


def compute_range(model_datas: Iterable[ModelData], function_data: FunctionData):
    #print("-----")
    #print(str(model_datas))
    #print(str(function_data))
    minimum, maximum = compute_boundaries(model_datas, function_data)
    return function_data.function_to_compute_ranges(maximum, minimum)


def validate_alleles(allele_in_model: str, alleles_in_data: Iterable[str], rsid: str, data_source: str):
    if allele_in_model not in alleles_in_data:
        print(str(allele_in_model) + " " + str(alleles_in_data))
        raise ConflictingAlleleBetweenDataAndModel(f'Conflicting alleles for {rsid} between model and {data_source}')


def validate_model_type(model_type: str, valid_model_types: Iterable[str]) -> str:
    small_letters_model_type = model_type.lower()
    if small_letters_model_type in valid_model_types:
        return small_letters_model_type
    raise NotImplementedError(f"Your model type: {small_letters_model_type} has not been implemented yet")


def compute_part_of_the_model(model_snips_and_coefficients: Dict[str, ModelData],
                              vcf_data: Union[Dict[str, SnpData], Dict[str, Tuple[str]]],
                              function_data: FunctionData) -> float:
    if not vcf_data:
        return function_data.neutral_number
    result = [None] * len(vcf_data)
    for i, (rsid, allele_data) in enumerate(vcf_data.items()):
        if isinstance(allele_data, SnpData):
            alleles = allele_data.genotype
        else:
            alleles = allele_data
        model_data = model_snips_and_coefficients[rsid]
        result[i] = function_data.function_to_apply_on_allele_pair(alleles, model_data.effect_allele,
                                                                   model_data.coeff_value)
    return function_data.function_to_apply_between_allele_pairs(result)


def rsids_data_from_data_source(rsids: Iterable[str], data_source: VcfAccessor, sample_name: str, is_imputed: bool = False) -> Dict[str, SnpData]:
    ret = {}
    for rsid in rsids:
        data = data_source.get_data_for_sample(sample_name, rsid, is_imputed)
        if data:
            ret[rsid] = data
    return ret


def get_allele_freq_in_poulation(snp_id: str, population_name: str,
                                 allele_frequency_accessor: AlleleFrequencyAccessor):
    logger.debug('Getting {} allele frequency from gnomad'.format(snp_id))
    return allele_frequency_accessor(snp_id, population_name)


def hardy_weinberg_diploid_frequencies(allele_frequencies: Dict[str, float]) -> Dict[Tuple[str], float]:
    assert len(allele_frequencies) == 2
    allele0 = list(allele_frequencies.keys())[0]
    allele1 = list(allele_frequencies.keys())[1]
    return {
        (allele0, allele0): allele_frequencies[allele0] ** 2,
        (allele1, allele1): allele_frequencies[allele1] ** 2,
        (allele0, allele1): allele_frequencies[allele1] * allele_frequencies[allele0] * 2
    }


def transform_boundaries_if_needed(model:PolygenicRiskScore, initial_boundaries:Dict[str, Dict[str, float]]):
    if model.model_type == 'beta':
        return initial_boundaries
    if model.model_type == 'odds_ratio':
       return change_boundaries(initial_boundaries)
    raise NotImplementedError(f'No known category boundaries transformation for model type {model.model_type}')


def change_boundaries(boundaries: Dict[str, Dict[str, float]]):
    ret_dict = {}
    for k, v in boundaries.items():
        ret_dict[k] = {}
        try:
            ret_dict[k]["from"] = math.log(v["from"])
        except ValueError:
            ret_dict[k]["from"] = np.NINF
        ret_dict[k]["to"] = math.log(v["to"])
    return ret_dict

if __name__ == '__main__':
    import glob, importlib, sys, os

    module_path = os.path.abspath(__file__).rsplit(os.path.sep, 6)[0]
    sys.path.insert(0, module_path)
    initial_paths = glob.glob('/data/test/large-data/vcf/mobigen_models_prepared/initial/*.vcf.gz')
    initial_accessor = VcfAccessor(initial_paths)
    imputed_paths = glob.glob('/data/test/large-data/vcf/mobigen_models_prepared/imputed/*.vcf.gz')
    imputed_accessor = VcfAccessor(imputed_paths)
    allele_accessor = AlleleFrequencyAccessor(
        allele_freq_json_path='/home/wojtek/PycharmProjects/mobigen/src/main/resources/allele_frequencies/gnomad.json')
    module = importlib.import_module('src/main/resources/test_models/testowy_nfe_model'.replace('/', '.'))
    model: PolygenicRiskScore = module.model
    data = Data(initial_accessor, imputed_accessor, allele_accessor, 'NA18500_9SR10670_D02.CEL', 'AF_nfe', model)
    print(data.compute_model())
    module2 = importlib.import_module('src/main/resources/test_models/testowy_covid_nfe_model'.replace('/', '.'))
    model2: PolygenicRiskScore = module2.model
    data2 = Data(initial_accessor, imputed_accessor, allele_accessor, 'NA18500_9SR10670_D02.CEL', 'AF_nfe', model2)
    print(data2.compute_model())
    module3 = importlib.import_module('src/main/resources/new_vitalleo_traits/covid_nfe_model'.replace('/', '.'))
    model3: PolygenicRiskScore = module3.model
    data3 = Data(initial_accessor, imputed_accessor, allele_accessor, 'NA18500_9SR10670_D02.CEL', 'AF_nfe', model3)
    print(data3.compute_model())
