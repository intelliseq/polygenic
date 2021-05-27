from src.main.python.lib.model.seqql import PolygenicRiskScore
from src.main.python.lib.model.seqql import ModelData
from src.main.python.lib.model.category import QuantitativeCategory

trait_was_prepared_for_population = "nfe"

model = PolygenicRiskScore(
    categories=[
        QuantitativeCategory(to=0.423, category_name='Low risk'),
        QuantitativeCategory(from_=0.423, to=1.01, category_name='Average risk'),
        QuantitativeCategory(from_=1.01, category_name='High risk')
    ],
    snips_and_coefficients={
        'rs35731912': ModelData(effect_allele='C', coeff_value=0.564000579501318),
        'rs657152': ModelData(effect_allele='C', coeff_value=1.32498349134725)
    },
    model_type='odds_ratio'
)
