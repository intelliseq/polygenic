import copy
from typing import Dict
from typing import Union
from typing import List
from typing import Sequence
from typing import Any


def vitalleo_model_output_single_gene(module) -> Dict[str, str]:
    return {
        'how_your_result_was_calculated': module.how_your_result_was_calculated
    }


def vitalleo_model_output_polygenic(module) -> Dict[str, Union[str, Sequence[str]]]:
    return {
        'genes': module.genes,
        'information_text': [module.how_your_result_was_calculated, module.science_behind_your_result]
    }


def vitalleo_snip(variants, patient_genotypes: Dict[str, List[str]],
                  report_numeric_gene_impact: bool = True):
    ret_variant_list = []
    for variant in variants:
        rsid = variant['rsid']
        patient_genotype = sorted(patient_genotypes[rsid])
        valid_patient_genotype = validate_snip_actionable(variant["genotypes"], patient_genotype, rsid)
        new_variant = copy.deepcopy(variant)
        new_variant["patient_genotype"] = valid_patient_genotype
        if report_numeric_gene_impact:
            new_variant['gene_impact'] = valid_patient_genotype.count(new_variant['risk_allele'])
        ret_variant_list.append(new_variant)
    return ret_variant_list


def validate_snip_actionable(genotypes: Dict[str, Any], patient_genotype: List[str], rsid: str):
    assert len(patient_genotype) in (1, 2)
    patient_genotype_set = set(patient_genotype)
    proper_genotypes = [set(genotype["genotype"]) for genotype in genotypes]
    if not patient_genotype_set in proper_genotypes:
        raise RuntimeError('Improper genotype: {} for rsid: {}. Proper genotypes are: {}'.format(
            patient_genotype, rsid, proper_genotypes))
    return patient_genotype if len(patient_genotype) == 2 else patient_genotype * 2


TEST_TYPE_TO_CALLABLE = {
    "Single SNP analysis": vitalleo_model_output_single_gene,
    "Polygenic Risk Score": vitalleo_model_output_polygenic
}


def vitalleo_model_output(category: str, module) -> Dict[str, Union[str, Sequence]]:
    ret_dict = {
        'trait': module.trait,
        'test_type': module.test_type,
        'about': module.about,
        'other_factors_statement': module.other_factors_statement,
        'other_factors': module.other_factors,
        'what_your_result_means': module.what_your_result_means[category],
        'result_statement': module.result_statement[category],
        'result': category,
        'what_you_can_do_statement': module.what_you_can_do_statement[category],
        'what_you_can_do': module.what_you_can_do[category],
        'references': module.references
    }
    ret_dict.update(TEST_TYPE_TO_CALLABLE[module.test_type](module))
    return ret_dict


def numeric_output(output: float, _) -> Dict[str, float]:
    return {
        'result': output
    }

def mobigen_output() -> Dict[str, float]:pass