import vcf
import sys
from polygenic.lib import polygenic
from polygenic.lib.seqql import QuantitativeCategory
from polygenic.lib.seqql import beta
from polygenic.lib.seqql import Choice
from polygenic.lib.seqql import Priority
from polygenic.lib.seqql import condition
from polygenic.lib.seqql import at_least
from polygenic.lib.seqql import Homozygote
from polygenic.lib.seqql import computation
from polygenic.lib.seqql import category
from polygenic.lib.seqql import always_true
from polygenic.lib.seqql import get_category
from polygenic.lib.seqql import odds_ratio

#################### TO DEFINE ####################
######### hair color
hair_color_model = Choice(
    {
        Priority(2):
            {
                # definicja warunku
                'when': at_least(
                    number_of_conditions_to_met=1,
                    conditions=[
                        condition(rs_id='rs1805005', genotype=Homozygote('T'), use_population_allele_frequency=False),
                        condition(rs_id='rs1805007', genotype=Homozygote('T'), use_population_allele_frequency=False),
                        condition(rs_id='rs1805009', genotype=Homozygote('C'), use_population_allele_frequency=False)
                    ],  # jeśli heterozygota, to Heterozygote('T/C')
                )
                        and
                        computation(
                            coefficients_and_snps_to_be_considered=[
                                beta(rs_id='rs12913832', effect_allele='A', beta_value=0.499),
                                beta(rs_id='rs16891982', effect_allele='C', beta_value=0.436),
                                beta(rs_id='rs12203592', effect_allele='T', beta_value=0.375),
                                beta(rs_id='rs1805007', effect_allele='T', beta_value=-0.356),
                                beta(rs_id='rs1426654', effect_allele='A', beta_value=-0.253),
                                beta(rs_id='rs12821256', effect_allele='T', beta_value=0.197),
                                beta(rs_id='rs6059655', effect_allele='A', beta_value=-0.182),
                                beta(rs_id='rs17184180', effect_allele='A', beta_value=-0.181),
                                beta(rs_id='rs72917317', effect_allele='T', beta_value=0.176),
                                beta(rs_id='rs80293268', effect_allele='C', beta_value=-0.171),
                                beta(rs_id='rs71443018', effect_allele='C', beta_value=-0.135)
                            ],
                            comparison_string='< 3',
                            function_to_apply='add'
                        ),

                # definicja akcji
                'action': category('red')  # zwrócenie dokładnej kategorii
            },
        Priority(1): {},  # opcjonalne akcje o niższym prioorytecie
        Priority(0):  # najniższy priorytet - defaultowe akcje,
        # jeśli nie zwróciliśmy żadnej kategorii dla wyższego priorytetu
            {
                'when': always_true,  # nie ma warunku, mamy akcje defaultowe:
                'action': get_category(
                    categories=[
                        QuantitativeCategory(to=0.5, category_name='light_blonde_hair'),
                        QuantitativeCategory(from_=0.5, to=2, category_name='blonde_hair'),
                        QuantitativeCategory(from_=2, to=3.5, category_name='light_brown_hair'),
                        QuantitativeCategory(from_=3.5, to=4.5, category_name='dark_brown_hair'),
                        QuantitativeCategory(from_=4.5, category_name='black_hair')
                    ],
                    coefficients_and_snps_to_be_considered=[
                        beta(rs_id='rs12913832', effect_allele='A', beta_value=0.499),
                        beta(rs_id='rs16891982', effect_allele='C', beta_value=0.436),
                        beta(rs_id='rs12203592', effect_allele='T', beta_value=0.375),
                        beta(rs_id='rs1805007', effect_allele='T', beta_value=-0.356),
                        beta(rs_id='rs1426654', effect_allele='A', beta_value=-0.253),
                        beta(rs_id='rs12821256', effect_allele='T', beta_value=0.197),
                        beta(rs_id='rs6059655', effect_allele='A', beta_value=-0.182),
                        beta(rs_id='rs17184180', effect_allele='A', beta_value=-0.181),
                        beta(rs_id='rs72917317', effect_allele='T', beta_value=0.176),
                        beta(rs_id='rs80293268', effect_allele='C', beta_value=-0.171),
                        beta(rs_id='rs71443018', effect_allele='C', beta_value=-0.135)
                    ],
                    function_to_apply='add'
                )
            }
    }
)
#### psoriasis
psoriasis_model = Choice(
    {
        Priority(0):
            {
                'when': always_true,  # nie ma dodatkowych warunków, jak przy kolore włosów
                'action': get_category(
                    categories=[
                        QuantitativeCategory(to=0.66, category_name='low_risk'),
                        QuantitativeCategory(from_=0.66, to=1.5, category_name='average_risk'),
                        QuantitativeCategory(from_=1.5, category_name='high_risk'),
                    ],
                    coefficients_and_snps_to_be_considered=[
                        odds_ratio(rs_id='rs35741374', effect_allele='A', odds_ratio_value=1.2),
                        odds_ratio(rs_id='rs1581803', effect_allele='A', odds_ratio_value=1.22),
                        odds_ratio(rs_id='rs643177', effect_allele='A', odds_ratio_value=1.27),
                        odds_ratio(rs_id='rs9487605', effect_allele='A', odds_ratio_value=1.27),
                        odds_ratio(rs_id='rs10515778', effect_allele='A', odds_ratio_value=1.29),
                        odds_ratio(rs_id='rs2546890', effect_allele='A', odds_ratio_value=1.39),
                        odds_ratio(rs_id='rs11135056', effect_allele='A', odds_ratio_value=1.45),
                        odds_ratio(rs_id='rs11575234', effect_allele='A', odds_ratio_value=1.47),
                        odds_ratio(rs_id='rs36207871', effect_allele='A', odds_ratio_value=1.47),
                        odds_ratio(rs_id='rs9481169', effect_allele='A', odds_ratio_value=1.58)
                    ],
                    function_to_apply='multiply'
                )
            }
    }
)
#################### END TO DEFINE ####################

if __name__ == '__main__':
    allele_frequency_in_population = {
        'rs9847240': {'A': 0.2, 'C': 0.3, 'G': 0.4, 'T': 0.1},
        'rs2025905': {'A': 0.25, 'C': 0.35, 'G': 0.3, 'T': 0.1},
        'rs1835873': {'A': 0.15, 'C': 0.3, 'G': 0.45, 'T': 0.1},
        'rs1805005': {'A': 0.8, 'C': 0.05, 'G': 0.05, 'T': 0.1}
    }  # TODO make sure this dictionary covers all snips and sum == 1
    # model = HairColourModel(allele_frequency_in_population)
    data = {}

    vcf_reader = vcf.Reader(sys.stdin)
    for record in vcf_reader:
        genotype = polygenic.get_genotype(record)
        data[record.ID] = genotype
        # print(data) 'rs79361800': 'G/C'
    # print([(x, x in data) for x in ['rs1805005', 'rs1805007', 'rs1805009']])
    # print([(x, x in data) for x in ['rs12913832', 'rs16891982', 'rs12203592', 'rs1805007', 'rs1426654', 'rs12821256',
    #                                 'rs6059655', 'rs17184180', 'rs72917317', 'rs80293268', 'rs71443018']])
    print(hair_color_model(data, allele_frequency_in_population))

    data = {'rs35741374': 'A/T', 'rs1581803': 'A/C', 'rs643177': 'A/C', 'rs9487605': 'A/C', 'rs10515778': 'A/C',
            'rs2546890': 'A/C', 'rs11135056': 'A/C', 'rs11575234': 'A/C', 'rs36207871': 'A/C', 'rs9481169': 'A/A'}
    print(psoriasis_model(data, allele_frequency_in_population))
