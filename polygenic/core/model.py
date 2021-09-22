import logging
import math
import yaml
from dotmap import DotMap

from polygenic.data.data_accessor import DataAccessor
from polygenic.core.utils import merge

logger = logging.getLogger('description_language.' + __name__)

class SeqqlOperator:

    def __init__(self, entries):
        self._entries = entries
        self.result = {}
        self._instantiate_subclass("model", Model)
        self._instantiate_subclass("description", Description)
        self._instantiate_subclass("categories", Categories)
        self._instantiate_subclass("diplotype_model", DiplotypeModel)
        self._instantiate_subclass("diplotypes", Diplotypes)
        self._instantiate_subclass("formula_model", FormulaModel)
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
        self.require("from")
        self.require("to")

    def compute(self, score: float):
        result = {"id": self.id, "match": False, "value": score}
        if score > self.get("from") and score < self.get("to"):
            result["match"] = True
        if self.has("scale_from") and self.has("scale_to"):
            result["value"] = self.get("scale_from") + (score - self.get("from")) / (self.get("to") - self.get("from")) * (self.get("scale_to") - self.get("scale_from"))
        return result
class DiplotypeModel(SeqqlOperator):
    def compute(self, data_accessor: DataAccessor):
        diplotypes_results = super(DiplotypeModel, self).compute(data_accessor)
        results = diplotypes_results["diplotypes"]
        results["genotypes"] = diplotypes_results["genotypes"]
        return results
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
            result["max"] += 2 * variant_result["effect_size"] if variant_result["effect_size"] > 0 else 0
            result["min"] += 2 * variant_result["effect_size"] if variant_result["effect_size"] < 0 else 0
            result["genotypes"][variant] = variant_result
            source = variant_result["genotype"]["source"]
            result["score"] += variant_result["score"]
            result[source + "_score"] += variant_result["score"]
            result[source + "_score_max"] += 2 * variant_result["effect_size"] if variant_result["effect_size"] > 0 else 0
            result[source + "_score_min"] += 2 * variant_result["effect_size"] if variant_result["effect_size"] < 0 else 0
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
        result["genotype"] = data_accessor.get_genotype_by_rsid(self.id)
        if self.has("diplotype"):
            result["diplotype_match"] = (sorted(self.get("diplotype").split('/')) == sorted(result["genotype"]["genotype"]))
        if self.has("effect_size") and self.has("effect_allele"):
            result["score"] = self.get("effect_size") * result["genotype"]["genotype"].count(self.get("effect_allele"))
            result["effect_size"] = self.get("effect_size")
        if self.has("symbol"):
            result["symbol"] = self.get("symbol")
        return result