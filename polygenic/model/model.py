import logging
import copy
import yaml
import math
from dotmap import DotMap

from polygenic.data.data_accessor import DataAccessor
from polygenic.model.utils import merge
from polygenic.error.polygenic_exception import PolygenicException

#tmp
import json

logger = logging.getLogger('description_language.' + __name__)

class SeqqlOperator:

    def __init__(self, entries):
        self._entries = entries
        self.__type__ = self.__class__.__name__
        self.result = {}
        self._instantiate_subclass("model", Model)
        self._instantiate_subclass("description", Description)
        self._instantiate_subclass("categories", Categories)
        self._instantiate_subclass("diplotype_model", DiplotypeModel)
        self._instantiate_subclass("diplotypes", Diplotypes)
        self._instantiate_subclass("formula_model", FormulaModel)
        self._instantiate_subclass("haplotype_model", HaplotypeModel)
        self._instantiate_subclass("haplotypes", Haplotypes)
        self._instantiate_subclass("score_model", ScoreModel)
        self._instantiate_subclass("variants", Variants)
        if type(self._entries) is dict:
            for entry in self._entries:
                if type(self._entries[entry]) is dict:
                    self._entries[entry] = SeqqlOperator(self._entries[entry])
        if type(self._entries) is list:
            for idx in range(len(self._entries)):
                if type(self._entries[idx]) is dict:
                    self._entries[idx] = SeqqlOperator(self._entries[idx])

    @classmethod
    def fromYaml(cls, path):
        seqql_yaml = {}
        with open(path, 'r') as stream:
            seqql_yaml = yaml.safe_load(stream)
        return cls(seqql_yaml)

    def _instantiate_subclass(self, key: str, cls):
        if (self._entries) is None:
            return
        for entry in self._entries:
            if entry == key:
                self._entries[key] = cls(self.get(key))

    def has(self, entry_name: str):
        if entry_name in self._entries:
            return True
        else:
            return False

    def get(self, entry_name: str):
        if self.has(entry_name):
            return(self._entries[entry_name])
        else:
            return None
    
    def set(self, entry_name: str, entry_value):
        self._entries[entry_name] = entry_value
        return self.has("entry_name")

    def get_entries(self):
        return self._entries
    
    def compute(self, data_accessor: DataAccessor):
        result = {}
        if type(self._entries) is not dict and type(self._entries) is not list:
            return self._entries
        if type(self._entries) is dict:
            for key in self._entries:
                if issubclass(self.get(key).__class__, SeqqlOperator):
                    merge(result, {key: self.get(key).compute(data_accessor)})
                else:
                    merge(result, {key: self.get(key)})
        if type(self._entries) is list:
            result = []
            for item in self._entries:
                #print(self.__type__ + " list ")
                if issubclass(item.__class__, SeqqlOperator):
                    result.append(item.compute(data_accessor))
        ### move genotypes to top
        for item in list(result):
            if type(result[item]) is dict and "genotypes" in result[item]:
                result["genotypes"] = result[item].pop("genotypes")
        self.result = result
        return result

    def require(self, key_name: str):
        if not self.has(key_name):
            raise RuntimeError(self.__class__.__name__ + " requires '" + key_name + "' field")

    def refine_results(self):
        """
        Selects category from fields marked with "_choice" based on results in model
        """
        refined_result = self.result
        category = "none"

        for item in self.result:
            if item.endswith("model") and "category" in self.result[item]:
                category = self.result[item]["category"]
        if "description" in self.result:
            description = self.result["description"]
            refined_description = {}
            for item in description:
                if item.endswith("_choice"):
                    if category is not None and category in description[item]:
                        refined_description[item[:-7]] = description[item][category]
                    else:
                        refined_description[item[:-7]] = None
                else:
                    refined_description[item] = description[item]
            refined_result["description"] = refined_description
        return refined_result

    def compute_qc(self, genotypes):
        qc = {"variant_count": 0,
              "variant_count_genotyping": 0,
                "variant_count_imputing": 0,
                "variant_count_af": 0,
                "variant_count_missing": 0,
                "variant_fraction_genotyping": 0,
                "variant_fraction_imputing": 0,
                "variant_fraction_af": 0,
                "variant_fraction_missing": 0,
              }
        for variant_id in genotypes:
            variant = genotypes[variant_id]
            qc["variant_count"] += 1
            if variant is not None and "source" in variant:
                if variant["source"] == "genotyping":
                    qc["variant_count_genotyping"] += 1
                elif variant["source"] == "imputing":
                    qc["variant_count_imputing"] += 1
                elif variant["source"] == "af":
                    qc["variant_count_af"] += 1
                elif variant["source"] == "missing":
                    qc["variant_count_missing"] += 1
                else:
                    raise PolygenicException("Unknown genotype source for " + variant_id)
            else:
                qc["variant_count_missing"] += 1
            if qc["variant_count"] > 0:
                qc["variant_fraction_genotyping"] = qc["variant_count_genotyping"] / qc["variant_count"]
                qc["variant_fraction_imputing"] = qc["variant_count_imputing"] / qc["variant_count"]
                qc["variant_fraction_af"] = qc["variant_count_af"] / qc["variant_count"]
                qc["variant_fraction_missing"] = qc["variant_count_missing"] / qc["variant_count"]
        return qc
class Description(SeqqlOperator):
    pass

class Categories(SeqqlOperator):
    def __init__(self, entries):
        super(Categories, self).__init__({})
        for key in entries:
            self._entries[key] = Category(key, entries[key])

class Category(SeqqlOperator):
    def __init__(self, key, entries):
        super(Category, self).__init__(entries)
        self.id = key

    def assign_category(self, category_name):
        result = {"id": self.id, "match": False, "category": self._entries}
        if (str(category_name) == str(self.id)):
            result["match"] = True
        return result

    def compute(self, score: float):
        result = {"id": self.id, "match": False, "value": score}
        if score >= self.get("from") and score <= self.get("to"):
            result["match"] = True
        if self.has("scale_from") and self.has("scale_to"):
            result["value"] = self.get("scale_from") + (score - self.get("from")) / (self.get("to") - self.get("from")) * (self.get("scale_to") - self.get("scale_from"))
        return result
    

class DiplotypeModel(SeqqlOperator):
    def compute(self, data_accessor: DataAccessor):
        #diplotypes_results = super(DiplotypeModel, self).compute(data_accessor)
        result = {}
        diplotype = self._entries["diplotypes"].compute(data_accessor)
        result["diplotype"] = diplotype["diplotype"]
        result["genotypes"] = diplotype["genotypes"]
        if self.has("categories"):
            for category_name in self.get("categories").get_entries():
                category = self.get("categories").get_entries()[category_name]
                category_result = category.assign_category(result["diplotype"])
                if category_result["match"]:
                    result["category"] = category_result["category"]
        result["qc"] = self.compute_qc(result["genotypes"])
        return result

class Diplotypes(SeqqlOperator):
    def __init__(self, entries):
        super(Diplotypes, self).__init__({})
        for diplotype in entries:
            self._entries[diplotype] = Diplotype(diplotype, entries[diplotype])

    def compute(self, data_accessor: DataAccessor):
        result = {"diplotype": None, "category": None}
        diplotypes_results = super(Diplotypes, self).compute(data_accessor)
        result["genotypes"] = diplotypes_results.pop("genotypes")
        for diplotype in diplotypes_results:
            if diplotypes_results[diplotype]["diplotype_match"]:
                result["diplotype"] = diplotype
                result["category"] = diplotype
        return result

class Diplotype(SeqqlOperator):
    def __init__(self, key, entries):
        super(Diplotype, self).__init__(entries)
        self.id = key

    def compute(self, data_accessor: DataAccessor):
        result = {}
        result["genotypes"] = {}
        result["diplotype_match"] = True
        variants_results = super(Diplotype, self).compute(data_accessor)["variants"]
        for variant in variants_results:
            variant_result = variants_results[variant]
            result["genotypes"][variant_result["genotype"]["rsid"]] = variant_result["genotype"]
            if not variant_result["diplotype_match"]:
                result["diplotype_match"] = False
        return result
class FormulaModel(SeqqlOperator):
    def __init__(self, entries):
        super(FormulaModel, self).__init__(entries)
        self.require("formula")

    def compute(self, data_accessor: DataAccessor):
        result = super(FormulaModel, self).compute(data_accessor)
        result["parameters"] = data_accessor.get_parameters()
        dotmap = DotMap(result)
        formula = self.get("formula").get_entries()
        for item in formula:
            expression = formula[item].replace("@", "dotmap.")
            value = eval(expression)
            dotmap[item] = value
            result[item] = value
        return result
class Model(SeqqlOperator):
    def __init__(self, entries):
        super(Model, self).__init__(entries)

class HaplotypeModel(SeqqlOperator):
    def compute(self, data_accessor: DataAccessor):
        if "haplotypes" not in self._entries:
            raise PolygenicException("HaplotypeModel requires genotypes field")
        genotypes = {}
        if "variants" in self._entries:
            genotypes = self._entries["variants"].compute(data_accessor)
        result = {}
        genotypes = self._entries["haplotypes"].compute_genotypes(data_accessor, genotypes)
        haplotypes = self._entries["haplotypes"].compute_haplotypes(genotypes)
        result["haplotypes"] = haplotypes
        result["genotypes"] = genotypes
        if self.has("categories"):
            for category_name in self.get("categories").get_entries():
                category = self.get("categories").get_entries()[category_name]
                category_result = category.compute(haplotypes["score"])
                if category_result["match"]:
                    result["category"] = category_name
                    result["value"] = category_result["value"]
        result["qc"] = self.compute_qc(result["genotypes"])
        return result

class Haplotypes(SeqqlOperator):
    def __init__(self, entries):
        super(Haplotypes, self).__init__({})
        for haplotype in entries:
            self._entries[haplotype] = Haplotype(haplotype, entries[haplotype])

    def compute_genotypes(self, data_accessor: DataAccessor, genotypes: dict):
        for entry in self._entries:
            genotypes = self._entries[entry].compute_genotypes(data_accessor, genotypes)
        return genotypes

    def compute_haplotypes(self, genotypes: dict):
        haplotypes = {}

        # compute best matching haplotype
        for haplotype_id in self._entries:
            haplotypes[haplotype_id] = self._entries[haplotype_id].compute_haplotype(genotypes)
            #print("COMPUTING 1: " + haplotype_id + " " + str(haplotypes[haplotype_id]['max_percent_match']))

        # sort and filter results by match
        sorted_haplotypes = sorted(haplotypes.items(), key = lambda x: x[1]['max_percent_match'], reverse = True)
        best_matching_haplotypes = {}
        FIRST_SORTED_ELEMENT = 0 # index of first element
        ID_IDX = 0
        HAPLOTYPE_IDX = 1 # haplotype is a tupple [ID, HAPLOTYPE]
        max_percent_match = sorted_haplotypes[FIRST_SORTED_ELEMENT][HAPLOTYPE_IDX]['max_percent_match']
        for haplotype in sorted_haplotypes:
            if haplotype[1]['max_percent_match'] >= (max_percent_match):
                best_matching_haplotypes[haplotype[ID_IDX]] = haplotype[HAPLOTYPE_IDX]

        # compute second haplotype
        matched_haplotypes = {}
        for first_haplotype_id in best_matching_haplotypes:
            computed_haplotypes = {}
            for second_haplotype_id in best_matching_haplotypes:
                leftover_genotypes = best_matching_haplotypes[first_haplotype_id]['leftover_genotypes'][0]
                computed_haplotypes[second_haplotype_id] = self._entries[second_haplotype_id].compute_haplotype(leftover_genotypes)
                #print("COMPUTING 2: " + first_haplotype_id + "_" + second_haplotype_id + " " + str(computed_haplotypes[second_haplotype_id]['max_percent_match']))
            sorted_haplotypes = sorted(computed_haplotypes.items(), key = lambda x: x[1]['max_percent_match'], reverse = True)
            max_percent_match = sorted_haplotypes[FIRST_SORTED_ELEMENT][HAPLOTYPE_IDX]['max_percent_match']
            sorted_haplotypes = sorted(computed_haplotypes.items(), key = lambda x: x[1]['min_percent_variants_missing'], reverse = False)
            for haplotype in sorted_haplotypes:
                if haplotype[1]['max_percent_match'] >= (max_percent_match):
                    id = first_haplotype_id + "_" + haplotype[0]
                    match = [best_matching_haplotypes[first_haplotype_id]['max_percent_match'],haplotype[1]['max_percent_match']]
                    missing = [best_matching_haplotypes[first_haplotype_id]['min_percent_variants_missing'],haplotype[1]['min_percent_variants_missing']]
                    first_haplotype = copy.deepcopy(best_matching_haplotypes[first_haplotype_id])
                    del first_haplotype['leftover_genotypes']
                    second_haplotype = copy.deepcopy(haplotype[1])
                    del second_haplotype['leftover_genotypes']
                    matched_haplotypes[id] = {
                        "id": id,
                        "match": match,
                        "missing": missing,
                        "match_sum": sum(match),
                        "missing_sum": sum(missing)
                    }
        sorted_haplotypes = sorted(matched_haplotypes.items(), key = lambda x: x[1]['match_sum'], reverse = True)
        match_sum = sorted_haplotypes[FIRST_SORTED_ELEMENT][HAPLOTYPE_IDX]['match_sum']
        match_filtered_haplotypes = {}
        for haplotype in sorted_haplotypes:
            if haplotype[1]['match_sum'] >= (match_sum):
                match_filtered_haplotypes[haplotype[0]] = haplotype[1] 
        sorted_haplotypes = sorted(match_filtered_haplotypes.items(), key = lambda x: x[1]['missing_sum'], reverse = False)
        missing_sum = sorted_haplotypes[FIRST_SORTED_ELEMENT][HAPLOTYPE_IDX]['missing_sum']
        missing_filtered_haplotypes = {}
        missing_unfiltered_haplotypes = {}
        for haplotype in sorted_haplotypes:
            if haplotype[1]['missing_sum'] <= (missing_sum):
                missing_filtered_haplotypes[haplotype[0]] = haplotype[1]
            else:
                missing_unfiltered_haplotypes[haplotype[0]] = haplotype[1]

        
        result = {}
        if max_percent_match >= 0.5:
            result["haplotype_id"] = list(missing_filtered_haplotypes.keys())[0]
        else:
            result["haplotype_id"] = None
        result["matching_haplotypes"] = [hap for hap in missing_filtered_haplotypes]
        result["unexcluded_haplotypes"] = [hap for hap in missing_unfiltered_haplotypes]
        result["id"] = result["haplotype_id"]
        result["haplotypes"] = matched_haplotypes
        return result

        
        result["id"] = [result["match"][0] + "/" + result["match"][1], result["match"][1] + "/" + result["match"][0]]
        result["score"] = score
        result["haplotypes"] = haplotypes
        return result

class Haplotype(SeqqlOperator):
    def __init__(self, key, entries):
        super(Haplotype, self).__init__({})
        self.id = key
        if entries:
            for key in entries:
                if key in ["af", "score"]:
                    self._entries[key] = entries[key]
                else:
                    self._entries[key] = Variant(key, entries[key])

    def compute_genotypes(self, data_accessor: DataAccessor, genotypes: dict):
        for entry in self._entries:
            if entry not in ["af", "score"] and entry not in genotypes:
                genotypes[entry] = self._entries[entry].compute(data_accessor)
        return genotypes

    def compute_haplotype(self, genotypes: dict):
        haplotype = {'id': self.id}
        required_matches = [0, 0] # total matches required
        matched_variants = [0, 0] # total variants matched
        missing_genotypes = [0, 0] # total missing variants
        missing_variants = [0, 0] # total missing haplotype specific variants
        percent_match = [0, 0] # percent of variants matched
        percent_genotypes_missing = [0, 0] # percent of genotypes missing
        percent_variants_missing = [0, 0] # percent of variants missing
        leftover_genotypes = [copy.deepcopy(genotypes), copy.deepcopy(genotypes)] # genotypes left to match in econd round

        # iterate over genotypes
        for genotype_id in genotypes:

            # increment required matches
            required_matches = [value + 1 for value in required_matches]

            # get genotype by id
            genotype = genotypes[genotype_id]

            # increment missing genotype counter
            if genotype["genotype"]["source"] == "missing":
                missing_genotypes = [value + 1 for value in missing_genotypes]
                if genotype_id in self._entries:
                    missing_variants = [value + 1 for value in missing_variants]

            # boolean allele match
            alleles = genotype["genotype"]["genotype"]
            if genotype_id in self._entries:
                alleles_match = [allele is not None and allele == genotype["effect_allele"] for allele in alleles]
            else:
                alleles_match = [allele is not None and allele != genotype["effect_allele"] for allele in alleles]
            if genotype["genotype"]["phased"]:
                matched_variants = [matched_variants[0] + alleles_match[0], matched_variants[1] + alleles_match[1]]
                leftover_genotypes[0][genotype_id]["genotype"]["genotype"] = [alleles[0], None]
                leftover_genotypes[1][genotype_id]["genotype"]["genotype"] = [alleles[1], None]
            else:
                if alleles_match[0]:
                    matched_variants = [matched_variants[0] + 1, matched_variants[1] + 1]
                    leftover_genotypes[0][genotype_id]["genotype"]["genotype"] = [alleles[1], None]
                    leftover_genotypes[1][genotype_id]["genotype"]["genotype"] = [alleles[1], None]
                elif alleles_match[1]:
                    matched_variants = [matched_variants[0] + 1, matched_variants[1] + 1]
                    leftover_genotypes[0][genotype_id]["genotype"]["genotype"] = [alleles[0], None]
                    leftover_genotypes[1][genotype_id]["genotype"]["genotype"] = [alleles[0], None]
                else:
                    leftover_genotypes[0][genotype_id]["genotype"]["genotype"] = [alleles[0], None]
                    leftover_genotypes[1][genotype_id]["genotype"]["genotype"] = [alleles[0], None]

        # compute percent match
        if (required_matches[0] - missing_genotypes[0]) > 0:
            percent_match[0] = matched_variants[0] / (required_matches[0] - missing_genotypes[0])
        if (required_matches[1] - missing_genotypes[1]) > 0:
            percent_match[1] = matched_variants[1] / (required_matches[1] - missing_genotypes[1])
        percent_match_order = sorted(range(len(percent_match)), key=lambda k: percent_match[k])
        percent_genotypes_missing[0] = missing_genotypes[0] / required_matches[0]
        percent_genotypes_missing[1] = missing_genotypes[1] / required_matches[1]
        percent_variants_missing[0] = missing_variants[0] / required_matches[0]
        percent_variants_missing[1] = missing_variants[1] / required_matches[1]
        haplotype["percent_match"] = percent_match
        haplotype["percent_genotypes_missing"] = percent_genotypes_missing
        haplotype["percent_variants_missing"] = percent_variants_missing
        haplotype["max_percent_match"] = percent_match[percent_match_order[1]]
        haplotype["min_percent_genotypes_missing"] = percent_genotypes_missing[percent_match_order[1]]
        haplotype["min_percent_variants_missing"] = percent_variants_missing[percent_match_order[1]]
        haplotype["required_matches"] = required_matches
        haplotype["matched_variants"] = matched_variants
        haplotype["missing_genotypes"] = missing_genotypes
        haplotype["missing_variants"] = missing_variants
        haplotype["leftover_genotypes"] = leftover_genotypes
        haplotype["af"] = self._entries["af"] if "af" in self._entries else 0
        haplotype["score"] = self._entries["score"] if "score" in self._entries else 0
        return haplotype


class ScoreModel(SeqqlOperator):
    def __init__(self, entries):
        super(ScoreModel, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        variants = self.get("variants").compute(data_accessor)
        result = {
            "value": 0,
            "constant": 0,
            "score": 0, 
            "max": 0,
            "min": 0,
            "genotypes": {}
        }
        for source in ["genotyping", "imputing", "af", "missing"]:
            result[source + "_score"] = 0
            result[source + "_score_max"] = 0
            result[source + "_score_min"] = 0
            result[source + "_alleles_count"] = 0
        for variant in variants:
            variant_result = variants[variant]
            effect_size = variant_result["effect_size"]
            result["max"] += 2 * effect_size if effect_size > 0 else 0
            result["min"] += 2 * effect_size if effect_size < 0 else 0
            result["genotypes"][variant] = variant_result
            source = variant_result["genotype"]["source"]
            result["score"] += variant_result["score"]
            result[source + "_score"] += variant_result["score"]
            result[source + "_score_max"] += 2 * effect_size if effect_size > 0 else 0
            result[source + "_score_min"] += 2 * effect_size if effect_size < 0 else 0
            result[source + "_alleles_count"] += 2            
        if self.has("constant"):
            result["score"] += self.get("constant")
            result["constant"] = self.get("constant")
        result["value"] = result["score"]
        if self.has("categories"):
            for category_name in self.get("categories").get_entries():
                category = self.get("categories").get_entries()[category_name]
                category_result = category.compute(result["score"])
                if category_result["match"]:
                    result["category"] = category_name
                    result["value"] = category_result["value"]
        result["qc"] = self.compute_qc(result["genotypes"])
        return result

class Variants(SeqqlOperator):
    def __init__(self, entries):
        super(Variants, self).__init__({})
        for key in entries:
            self._entries[key] = Variant(key, entries[key])
    def compute(self, data_accessor:DataAccessor):
        return super(Variants, self).compute(data_accessor)

class Variant(SeqqlOperator):
    def __init__(self, key, entries):
        super(Variant, self).__init__(entries)
        self.id = key

    def compute(self, data_accessor: DataAccessor):
        result = {}
        if not self._entries:
            result["genotype"] = {'rsid': self.id, 'genotype': [None, None], 'phased': None, 'source': 'invalidmodelentry', 'ref': None}
            result["effect_allele"] = None
            return result
        result["genotype"] = data_accessor.get_genotype_by_rsid(self.id)
        if self.has("diplotype"):
            if result["genotype"]["genotype"][0] is None:
                result["diplotype_match"] = False
            else:
                result["diplotype_match"] = (sorted(self.get("diplotype").split('/')) == sorted(result["genotype"]["genotype"]))
            result["effect_allele"] = None
        if self.has("effect_size") and self.has("effect_allele"):

            ### CHECK effect_size
            try: effect_size = float(self.get("effect_size"))
            except: raise PolygenicException("bad or missing effect_size for variant: " + str(self._entries))
            self.set("effect_size", effect_size)

            result["score"] = self.get("effect_size") * result["genotype"]["genotype"].count(self.get("effect_allele"))
            result["effect_size"] = self.get("effect_size")
        if self.has("symbol"):
            result["symbol"] = self.get("symbol")
        if self.has("effect_allele"):
            result["effect_allele"] = self.get("effect_allele")
        elif self.has("alt"):
            result["effect_allele"] = self.get("alt")
        else:
            result["effect_allele"] = None
        return result