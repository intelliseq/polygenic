import logging
import math
import yaml
from dotmap import DotMap

from polygenic.lib.data_access.data_accessor import DataAccessor
from polygenic.lib.utils import merge

logger = logging.getLogger('description_language.' + __name__)

class SeqqlOperator:

    def __init__(self, entries):
        self._entries = entries
        self._define_parameter("model", Model)
        self._define_parameter("score", Score)
        self._define_parameter("exp", Exp)
        self._define_parameter("label", Label)
        self._define_parameter("eval", Eval)
        self._define_parameter("items", Items)
        self._define_parameter("diplotype_model", DiplotypeModel)
        self._define_parameter("diplotype_variants", DiplotypeVariants)
        self._define_parameter("diplotypes", Diplotypes)
        self._define_parameter("variants", Variants)
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
    
    def compute(self, data_accessor: DataAccessor):
        result = {}
        if type(self._entries) is dict:
            for key in self._entries:
                if issubclass(self.get(key).__class__, SeqqlOperator):
                    merge(result, self.get(key).compute(data_accessor))
        if type(self._entries) is list:
            for item in self._entries:
                if issubclass(item.__class__, SeqqlOperator):
                    merge(result, item.compute(data_accessor))
        return result

    def require(self, key_name: str):
        if not self.has(key_name):
            raise RuntimeError(self.__class__.__name__ + " requires '" + key_name + "' component")

    def _define_parameter(self, key: str, cls):
        for entry in self._entries:
            if entry == key:
                self._entries[key] = cls(self.get(key))

class Eval(SeqqlOperator):
    def __init__(self, entries):
        super(Eval, self).__init__(entries)
        self.require("formula")

    def compute(self, data_accessor: DataAccessor):
        computation_result = super(Eval, self).compute(data_accessor)
        dotmap = DotMap(computation_result)
        formula = self.get("formula")
        formula = formula.replace("@", "dotmap.")
        computation_result["formula"] = formula
        computation_result["value"] = eval(formula)
        return computation_result

class Exp(SeqqlOperator):
    def __init__(self, entries):
        super(Exp, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        result = super(Exp, self).compute(data_accessor)
        if "value" in result:
            result["value"] = math.exp(result["value"])
        return result

class DiplotypeModel(SeqqlOperator):
    def __init__(self, entries):
        super(DiplotypeModel, self).__init__(entries)
class Diplotypes(SeqqlOperator):
    def __init__(self, entries):
        super(Diplotypes, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        #result = super(Diplotypes, self).compute(data_accessor)
        result = {}
        for diplotype in self._entries:
            Diplotype
            #diplotype_result = diplotype.compute(data_accessor)
            #result[diplotype["name"]] = diplotype_result
        return result
class DiplotypeVariants(SeqqlOperator):
    def __init__(self, entries):
        super(DiplotypeVariants, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        results = {
            "value": False,
            "genotyping_alleles_count": 0,
            "imputing_alleles_count": 0,
            "af_alleles_count": 0, 
            "missing_alleles_count": 0,
            "genotypes": []
        }
        for variant in self._entries:
            genotype = data_accessor.get_genotype_by_rsid(variant.get("rsid"))
            if (genotype["source"]) == "missing":
                results["value"] = False
            else:
                if (genotype["genotype"].sort()[0] == variant.get("alleles")[0]) and (genotype["genotype"].sort()[1] == variant.get("alleles")[1]):
                    results["value"] = True
        return results
class Items(SeqqlOperator):
    def __init__(self, entries):
        super(Items, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        result = super(Items, self).compute(data_accessor)
        return result
class Model(SeqqlOperator):
    def __init__(self, entries):
        super(Model, self).__init__(entries)
class Label(SeqqlOperator):
    def __init__(self, entries):
        super(Label, self).__init__(entries)
        self.require("name")

    def compute(self, data_accessor: DataAccessor):
        label_results = {}
        genotypes = []
        computation_results = super(Label, self).compute(data_accessor)
        if "genotypes" in computation_results:
            for genotype in computation_results["genotypes"]:
                genotype["label"] = self.get("name")
                genotypes.append(genotype)
            label_results["genotypes"] = genotypes
            del computation_results["genotypes"]
        computation_results["label"] = self.get("name")
        if ("label" in computation_results) and ("value" in computation_results):
            label_results[computation_results["label"]] = computation_results
            label_results["value"] = computation_results["value"]
        return label_results
class Score(SeqqlOperator):
    def __init__(self, entries):
        super(Score, self).__init__(entries)


    def compute(self, data_accessor: DataAccessor):
        computation_results = {"result": 0}
        if self.has("label"): computation_results.update(self.get("label").compute(data_accessor))
        if self.has("variants"): computation_results.update(self.get("variants").compute(data_accessor))
        return computation_results
        #return {"z": data_accessor.get_genotype_by_rsid("rs10012")}

class Variants(SeqqlOperator):
    def __init__(self, entries):
        super(Variants, self).__init__(entries)

    def compute(self, data_accessor: DataAccessor):
        computation_results = {
            "value": 0,
            "constant": 0,
            "score": 0, 
            "cap": 0,
            "genotyping_score": 0, 
            "imputing_score": 0, 
            "af_score": 0,
            "missing_score": 0,
            "genotyping_score_cap": 0,
            "imputing_score_cap": 0,
            "af_score_cap": 0, 
            "missing_score_cap": 0,
            "genotyping_alleles_count": 0,
            "imputing_alleles_count": 0,
            "af_alleles_count": 0, 
            "missing_alleles_count": 0,
            "genotypes": []
        }
        for variant in self._entries:
            genotype = data_accessor.get_genotype_by_rsid(variant.get("rsid"))
            computation_results["cap"] = computation_results["cap"] + 2 * variant.get("effect_size")
            if len(genotype.keys()) < 2:
                genotype["genotype"] = [None, None]
                genotype["source"] = "missing"
                computation_results["genotypes"].append(genotype)
                computation_results["missing_score_cap"] += 2 * variant.get("effect_size")
                computation_results["missing_alleles_count"] += 2
            else:
                source = genotype["source"]
                computation_results["genotypes"].append(genotype)
                for allele in genotype["genotype"]:
                    if allele == variant.get("effect_allele"):
                        computation_results["score"] += variant.get("effect_size")
                        computation_results[source + "_score"] += variant.get("effect_size")
                computation_results[source + "_score_cap"] += 2 * variant.get("effect_size")
                computation_results[source + "_alleles_count"] += 2
            computation_results["value"] = computation_results["score"]
            if self.has("constant"):
                computation_results["value"] += self.get("constant")
                computation_results["constant"] = self.get("constant")
        return computation_results