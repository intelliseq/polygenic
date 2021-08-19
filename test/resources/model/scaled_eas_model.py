from polygenic.seqql.score import PolygenicRiskScore
from polygenic.seqql.score import ModelData
from polygenic.seqql.category import QuantitativeCategory

trait_was_prepared_for_population = 'eas'

model = PolygenicRiskScore(
    categories=[
        QuantitativeCategory(from_=1.371624087, to=2.581880425, category_name='High risk', scale_from = 0, scale_to = 1),
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

description = {
  'about': '',
  'genes': [],
  'result_statement_choice': {
    'Average risk': 'Avg',
    'Potential risk': 'Pot',
    'High risk': 'Hig',
    'Low risk': 'Low'
  },
  'science_behind_the_test': '',
  'test_type': 'Polygenic Risk Score',
  'trait': 'Breast cancer',
  'trait_authors': [
    'taken from the PGS catalog'
  ],
  'trait_copyright': 'Intelliseq all rights reserved',
  'trait_explained': None,
  'trait_heritability': None,
  'trait_pgs_id': 'PGS000001',
  'trait_pmids': [
    '25855707'
  ],
  'trait_snp_heritability': None,
  'trait_title': 'Breast_Cancer',
  'trait_version': '1.0',
  'what_you_can_do_choice': {
    'Average risk': '',
    'High risk': '',
    'Low risk': ''
  },
  'what_your_result_means_choice': {
    'Average risk': '',
    'High risk': '',
    'Low risk': ''
  }
}