import logging
import math
import yaml
from dotmap import DotMap

from polygenic.data.data_accessor import DataAccessor
from polygenic.model.utils import merge
from polygenic.error import polygenic_exception

#tmp
import json

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
        if self._entries and not isinstance(self._entries, list):
            self.require("from")
            self.require("to")

    def compute(self, score: float):
        result = {"id": self.id, "match": False, "value": score}
        if score >= self.get("from") and score <= self.get("to"):
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

class HaplotypeModel(SeqqlOperator):
    def compute(self, data_accessor: DataAccessor):
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
        matching_haplotypes = {}
        result = {}
        for entry in self._entries:
            haplotypes[entry] = self._entries[entry].compute_haplotype(genotypes)
        for haplotype_id in haplotypes:
            haplotype = haplotypes[haplotype_id]
            if (sum(haplotype["alleles"]) > 0):
                matching_haplotypes[haplotype_id] = {
                    "id": haplotype_id,
                    "allele_count": sum(haplotype["alleles"]),
                    "af": haplotype["af"],
                    "alleles": haplotype["alleles"],
                    "score": haplotype["score"]
                }
        res = sorted(matching_haplotypes.items(), key = lambda x: x[1]['af'])
        allele_count_sum = 0
        result["match"] = []
        score = 0
        for item in res:
            haplotype_id = item[0]
            haplotype = item[1]
            for i in range(int(haplotype["allele_count"])):
                score += haplotype["score"]
                result["match"].append(haplotype["id"])
                allele_count_sum += 1
                if allele_count_sum >= 2: break
            if allele_count_sum >= 2: break
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
        haplotype_tupple = [True, True]
        possible_haplotype_allele_count = 2
        required_matches = [0, 0]
        matched_variants = [0, 0]
        missing_genotypes = 0
        missing_variants = 0
        for genotype_id in genotypes:
            required_matches = [required_matches[0] + 1, required_matches[1] + 1]
            genotype = genotypes[genotype_id]
            if genotype["genotype"]["source"] == "missing":
                missing_genotypes += 1
                if genotype_id in self._entries:
                    missing_variants += 1
            alleles = [allele == genotype["effect_allele"] for allele in genotype["genotype"]["genotype"]]
            if genotype_id not in self._entries:
                alleles = [not allele for allele in alleles]
            phased = genotype["genotype"]["phased"]
            allele = alleles[0] or alleles[1]
            possible_haplotype_allele_count = min(possible_haplotype_allele_count, sum(alleles))
            if phased:
                matched_variants = [matched_variants[0] + alleles[0], matched_variants[1] + alleles[1]]
                haplotype_tupple[0] = haplotype_tupple[0] and alleles[0]
                haplotype_tupple[1] = haplotype_tupple[1] and alleles[1]
            else:
                matched_variants = [matched_variants[0] + allele, matched_variants[1] + allele]
                haplotype_tupple[0] = haplotype_tupple[0] and allele
                haplotype_tupple[1] = haplotype_tupple[1] and allele
        haplotype["possible_haplotype_allele_count"] = possible_haplotype_allele_count
        haplotype["required_matches"] = required_matches
        haplotype["matched_variants"] = matched_variants
        haplotype["missing_genotypes"] = missing_genotypes
        haplotype["missing_variants"] = missing_variants
        haplotype["alleles"] = haplotype_tupple
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
            result["effect_allele"] = self.get("alt")
            return result
        result["genotype"] = data_accessor.get_genotype_by_rsid(self.id)
        if self.has("diplotype"):
            result["diplotype_match"] = (sorted(self.get("diplotype").split('/')) == sorted(result["genotype"]["genotype"]))
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