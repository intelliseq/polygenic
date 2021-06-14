import functools
import abc
import operator
import inspect
import logging
import re
import numpy as np
from typing import Callable
from typing import Iterable
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from polygenic.lib.allele_frequency_accessor import AlleleFrequencyAccessor

logger = logging.getLogger('description_language.' + __name__)

FUNCTION_DICT = {
    'add': sum,
    'multiply': functools.partial(functools.reduce, operator.mul)
}

FUNCTION_NAME_RE = re.compile('function (.*)\.<locals>')


class GenotypeGetter(object):

    def __init__(self, rsids: List) -> None:
        super().__init__()
        self.rsids = rsids

    def __call__(self, data):
        return {rsid: data[rsid] for rsid in self.rsids}


class Zygosity(abc.ABC):
    def __init__(self, allele):
        super().__init__()
        self.allele = allele

    @abc.abstractmethod
    def __eq__(self, other_genotype):
        if not (isinstance(other_genotype, list) or isinstance(other_genotype, tuple)):
            return False


class Homozygote(Zygosity):
    def __init__(self, allele: str):
        assert isinstance(allele, str)
        super().__init__(allele)

    def __eq__(self, other_genotype: Union[List, Tuple]):
        super().__eq__(other_genotype)
        other_genotype_set = set(other_genotype)
        return len(other_genotype_set) == 1 and self.allele in other_genotype_set


class Heterozygote(Zygosity):
    def __init__(self, alleles):
        assert isinstance(alleles, list) or isinstance(alleles, tuple) or isinstance(alleles, str)
        if isinstance(alleles, str):
            alleles = alleles.split('/')
        assert len(alleles) == 2
        super().__init__(alleles)

    def __eq__(self, other_genotype):
        super().__eq__(other_genotype)
        other_genotype_set = set(other_genotype)
        return len(other_genotype_set) == 2 and set(self.allele) == other_genotype_set


def at_least(number_of_conditions_to_met: int, conditions: List[Callable]):
    def f(data: Dict[str, Union[List, Tuple]], population_name: str,
          allele_frequency_accessor: AlleleFrequencyAccessor):
        conditions_met = 0
        for condition in conditions:
            if condition(data, population_name, allele_frequency_accessor):
                conditions_met += 1
            if conditions_met >= number_of_conditions_to_met:
                return True
        return False

    return f


def condition(rs_id, genotype, allele_frequency_accessor: AlleleFrequencyAccessor,
              use_population_allele_frequency=False) -> bool:
    def f(data: Dict[str, str], population_name: str):
        # data - słownik postaci rsid:'A/C'
        alleles = get_allele_frequency(rs_id, data, use_population_allele_frequency, population_name,
                                       allele_frequency_accessor)
        return genotype == alleles

    return f


def get_allele_frequency(rs_id: str, data: Dict[str, str], use_population_allele_frequency: bool, population_name: str,
                         allele_frequency_accessor: AlleleFrequencyAccessor):
    if use_population_allele_frequency:
        return get_alleles_when_missing_consider_population_freqs(rs_id, data, population_name,
                                                                  allele_frequency_accessor)
    else:
        return data[rs_id]


def get_value(coefficients_and_snps_to_be_considered: List[Callable], function_to_apply: str,
              function_dict: Dict[str, Callable] = FUNCTION_DICT):
    def f(data, population: str, allele_frequency_accessor: AlleleFrequencyAccessor):
        s = function_dict[function_to_apply](
            snp_function(data, population, allele_frequency_accessor) for snp_function in
            coefficients_and_snps_to_be_considered)
        return s

    return f


def get_category_and_value(categories: List, coefficients_and_snps_to_be_considered: List[Callable],
                           function_to_apply: str,
                           function_dict: Dict[str, Callable] = FUNCTION_DICT) -> Dict[str, Union[str, float]]:
    """
    Factory function that returns a function to compute category and numerical output for given coefficients and genotype
    :param categories: list of possible categories
    :param coefficients_and_snps_to_be_considered: list of callables to be applied to the genotype
    :param function_to_apply: string representation of function to be applied to the coefficients
    :param function_dict: translation between string representation of a function and the corresponding python callable
    :return: function to compute category and numerical output for given coefficients and genotype
    """

    def f(data, population: str, allele_frequency_accessor: AlleleFrequencyAccessor):
        s = function_dict[function_to_apply](
            snp_function(data, population, allele_frequency_accessor) for snp_function in
            coefficients_and_snps_to_be_considered)
        for category in categories:
            if category.from_ <= s <= category.to:
                returned_category_name = category.name
                break
        return {'categorical_output': returned_category_name, 'numerical_output': s}

    return f


def computation(coefficients_and_snps_to_be_considered: List[Callable], comparison_string: str, function_to_apply: str,
                function_dict: Dict[str, Callable] = FUNCTION_DICT):
    def f(data: Dict[str, str], population_name: str):
        return eval('{}{}'.format(function_dict[function_to_apply](
            snp_function(data, population_name) for snp_function in coefficients_and_snps_to_be_considered),
            comparison_string))

    return f


def odds_ratio(rs_id: str, effect_allele: str, odds_ratio_value: float) -> Callable:
    def func(data, population_name: str, allele_frequency_accessor: AlleleFrequencyAccessor):
        alleles = get_alleles_when_missing_consider_population_freqs(rs_id, data, population_name,
                                                                     allele_frequency_accessor)
        if len(alleles) < 1:
            logger.error('Improper genotype for {}: {}'.format(rs_id, alleles))
            return 1
        if alleles[0] == alleles[1] == effect_allele:
            return odds_ratio_value
        elif effect_allele in alleles:
            return 1
        else:
            return 1 / odds_ratio_value

    return func


class Priority(int):
    """
    Human-friendly representation of priority of a condition
    """

    def __new__(cls, priority: int):
        return super().__new__(cls, priority)

    def __gt__(self, other):
        if isinstance(other, self.__class__):
            return super().__gt__(self, other)
        else:
            raise TypeError('unorderable types: {} and {}'.format(type(self), type(other)))


class Choice(object):
    '''
    Model type which checks some conditions according to their priorities. When a condition is met, the part of the model for this condition is evaluated
    '''

    def __init__(self, priority_conditions_action_dict: Dict[Priority, Dict[str, Callable]]):
        """
        Initializes Choice instance
        :param priority_conditions_action_dict: dictionary of type {priority: {'when': function_that_tells_when_the_condition_is_evaluated, 'action': function_that_tells_how_to_evaluate_a_given_condition}}
        """
        super().__init__()
        self.priority_conditions_action_dict = priority_conditions_action_dict

    def __call__(self, data: Dict[str, str], allele_frequency_accessor: AlleleFrequencyAccessor, population):
        # data - a dictionary rsid:'A/C'
        for priority, condition_action_dict in sorted(self.priority_conditions_action_dict.items(), reverse=True):
            if condition_action_dict['when'](data, population, allele_frequency_accessor):
                return condition_action_dict['action'](data, population, allele_frequency_accessor)
        raise NoCategoryError("None od possible categories matched analysed data")

    def get_rsids_used_by_model(self) -> List[str]:
        # does not work for al_least function!
        callables = []
        for priority_value in self.priority_conditions_action_dict.values():
            if priority_value:
                # print('prio_val', priority_value['when'], inspect.getclosurevars(priority_value['when']).nonlocals)
                action_callables = self.__get_callable_from_priority_conditions_action_dict_val(priority_value,
                                                                                                'action')
                if action_callables:
                    callables.extend(action_callables)
                when_callables = self.__get_callable_from_priority_conditions_action_dict_val(priority_value, 'when')
                if when_callables:
                    callables.extend(when_callables)
        # print(callables)
        return [self.__get_rsid_from_callable(callable) for callable in callables]

    def category_boundaries(self, priority=Priority(0)) -> Dict[str, Dict[str, str]]:
        if priority != 0:
            raise NotImplementedError("Implmented only for priority = 0")
        categories = inspect.getclosurevars(self.priority_conditions_action_dict[priority]['action']).nonlocals[
            'categories']
        return {category.name: {'from': category.from_, 'to': category.to} for category in categories}

    def model_type(self) -> str:
        names_of_functions_used_to_compute_final_result = set()
        for _, v in self.priority_conditions_action_dict.items():
            function_name = FUNCTION_NAME_RE.search(str(v['action'])).group(1)
            if not function_name in ('get_value', 'get_category_and_value'):
                raise NotImplementedError('Not implemented for "{}" function'.format(function_name))
            for x in inspect.getclosurevars(v['action']).nonlocals['coefficients_and_snps_to_be_considered']:
                names_of_functions_used_to_compute_final_result.add(FUNCTION_NAME_RE.search(str(x)).group(1))
        assert len(names_of_functions_used_to_compute_final_result) == 1
        function_name = names_of_functions_used_to_compute_final_result.pop()
        return function_name

    @staticmethod
    def __get_callable_from_priority_conditions_action_dict_val(priority_conditions_action_dict_val: Dict, key: str) -> \
            Union[List[Callable], None]:
        try:
            return inspect.getclosurevars(
                priority_conditions_action_dict_val[key]).nonlocals['coefficients_and_snps_to_be_considered']
        except KeyError:
            pass

    @staticmethod
    def __get_rsid_from_callable(callale_with_rsid_in_closurevars: Callable) -> str:
        '''
        Method to extract rs_id from functions such as 'beta', 'odds_ratio', 'condition'
        :param callale_with_rsid_in_closurevars:
        :return:
        '''
        return inspect.getclosurevars(callale_with_rsid_in_closurevars).nonlocals['rs_id']


class NoCategoryError(RuntimeError):
    pass


def always_true(*args, **kwargs) -> bool:
    """
    Always returns true
    :param args:
    :param kwargs:
    :return: True
    """
    return True


def category(category_string: str) -> Callable:
    def f(_, __):
        return category_string

    return f


class QuantitativeCategory(object):
    def __init__(self, category_name, from_=-np.inf, to=np.inf):
        super().__init__()
        self.from_ = from_
        self.to = to
        self.name = category_name

    def intersects_with(self, other) -> None:
        if self.from_ == other.from_ or self.to == other.to:
            return True
        elif self.from_ < other.from_ < self.to:
            return True
        elif other.from_ < self.from_ < other.to:
            return True
        return False


class InvalidCategoriesError(RuntimeError):
    pass


class NoCategoryError(RuntimeError):
    """
    Raised  when no category could be returned.
    """
    pass


def get_alleles_when_missing_consider_population_freqs(rs_id: str, data: Dict[str, Union[List, Tuple]],
                                                       population_name: str,
                                                       allele_frequency_accessor: AlleleFrequencyAccessor):
    try:
        ret_val = data[rs_id]
        return ret_val
    except KeyError:
        logger.debug('{} was not found in input data'.format(rs_id))
        return max(hardy_weinberg_diploid_frequencies(
            get_allele_freq_in_poulation(rs_id, population_name, allele_frequency_accessor)).items(),
                   key=lambda x: x[1])[0]


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


def beta(rs_id: str, effect_allele: str, beta_value: float) -> Callable:
    def func(data, population_name: str, allele_frequency_accessor: AlleleFrequencyAccessor):
        alleles = get_alleles_when_missing_consider_population_freqs(rs_id, data, population_name,
                                                                     allele_frequency_accessor)
        result = 0
        for allele in alleles:
            if allele == effect_allele:
                result += beta_value
        return result

    return func


def validate_categories(categories: Iterable[QuantitativeCategory]) -> None:
    validate_categories_names(categories)
    validate_categories_ranges(categories)


def validate_categories_names(categories: Iterable[QuantitativeCategory]):
    sorted_according_to_names = sorted(categories, key=lambda x: x.name)
    for i in range(len(sorted_according_to_names) - 1):
        if categories_have_same_names(sorted_according_to_names[i], sorted_according_to_names[i + 1]):
            raise InvalidCategoriesError(
                'Categories {} and {} have the same name'.format(sorted_according_to_names[i],
                                                                 sorted_according_to_names[i + 1]))


def categories_have_same_names(category1: QuantitativeCategory, category2: QuantitativeCategory):
    return category1.name == category2.name


def validate_categories_ranges(categories: Iterable[QuantitativeCategory]):
    for i in range(len(categories) - 1):
        for j in range(1, len(categories)):
            if categories[i].intersects_with(categories[j]):
                raise InvalidCategoriesError('Categories {} and {} intersect'.format(categories[i], categories[j]))


if __name__ == '__main__':
    print(Homozygote('A') == ['A', 'A'])
    het = Heterozygote('A/C')
    print(het.__eq__(['A', 'C']))
    print(het == ['C', 'A'])
    print(het == ['A', 'A'])
    print(het == ['T', 'T'])
    print(het == ['C'])
