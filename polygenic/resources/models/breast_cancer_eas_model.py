from polygenic.lib.model.seqql import PolygenicRiskScore
from polygenic.lib.model.seqql import ModelData
from polygenic.lib.model.category import QuantitativeCategory

trait_was_prepared_for_population = "eas"

model = PolygenicRiskScore(
    categories=[
        QuantitativeCategory(from_=1.371624087, to=2.581880425, category_name='High risk'),
        QuantitativeCategory(from_=1.169616034, to=1.371624087, category_name='Potential risk'),
        QuantitativeCategory(from_=-0.346748358, to=1.169616034, category_name='Average risk'),
	QuantitativeCategory(from_=-1.657132197, to=-0.346748358, category_name='Low risk')
    ],
    snips_and_coefficients={
	'rs10012': ModelData(effect_allele='G', coeff_value=0.369215857410143),
	'rs1014971': ModelData(effect_allele='T', coeff_value=0.075546961392531),
	'rs10936599': ModelData(effect_allele='C', coeff_value=0.086359830674748),
	'rs11892031': ModelData(effect_allele='C', coeff_value=-0.552841968657781),
	'rs1495741': ModelData(effect_allele='A', coeff_value=0.05307844348342),
	'rs17674580': ModelData(effect_allele='C', coeff_value=0.187520720836463),
	'rs2294008': ModelData(effect_allele='T', coeff_value=0.08278537031645),
	'rs798766': ModelData(effect_allele='T', coeff_value=0.093421685162235),
	'rs9642880': ModelData(effect_allele='G', coeff_value=0.093421685162235)
    },
    model_type='beta'
)
