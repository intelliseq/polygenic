from polygenic.seqql import QuantitativeCategory
from polygenic.seqql import beta
from polygenic.seqql import Choice
from polygenic.seqql import Priority
from polygenic.seqql import always_true
from polygenic.seqql import get_category

######### type_2_diabetes
#na podstawie: doi:10.1038/ng.3943
#dane tylko dla populacji EUR, w sumie 4258 "cases" i 18519 "controls" (SuppTable1)
#dane o rs: z tabeli SuppTable5 dla populacji EUR; 
#wybralam warianty p<0.5e-08 i oddalone od innych istotnych wariantow o wiecej niż 500 tys.pz (zostawialam wariant o wiekszym efekcie);
#w sumie zostało ich 41
#efekty byly podane jako wspolczynniki beta z regresji logistycznej
#odziedziczalnosc: 20-80%, zmiennosc wyjasniona: nie podano

#model na podstawie 100000 symulacji (plik symul.py), wyniki w pliku dist_2td.txt i stst_t2d.txt (tu z pieciu powtorzen symulacji)
#czestosci dla EA z Ensembl (EUR, 1000Genomes)

type_2_diabetes_model = Choice(
    {
        Priority(0):  # najniższy priorytet - defaultowe akcje,
        # jeśli nie zwróciliśmy żadnej kategorii dla wyższego priorytetu
            {
                'when': always_true,  # nie ma warunku, mamy akcje defaultowe:
                'action': get_category(
                    categories=[
                        QuantitativeCategory(to=-0.25, category_name='low_risk'),
                        QuantitativeCategory(from_=-0.25, to=0.77, category_name='average_risk'),
                        QuantitativeCategory(from_=0.77, category_name='high_risk'),
                        
                    ],
                    coefficients_and_snps_to_be_considered=[
                        beta(rs_id='rs10811661', effect_allele='C', beta_value=-0.156),
                        beta(rs_id='rs10830963', effect_allele='G', beta_value=0.094),
                        beta(rs_id='rs10842994', effect_allele='T', beta_value=-0.081),
                        beta(rs_id='rs11063069', effect_allele='G', beta_value=0.074),
                        beta(rs_id='rs1111875', effect_allele='T', beta_value=-0.107),
                        beta(rs_id='rs11257655', effect_allele='T', beta_value=0.077),
                        beta(rs_id='rs11708067', effect_allele='G', beta_value=-0.104),
                        beta(rs_id='rs12571751', effect_allele='G', beta_value=-0.066),
                        beta(rs_id='rs12970134', effect_allele='A', beta_value=0.071),
                        beta(rs_id='rs13266634', effect_allele='T', beta_value=-0.111),
                        beta(rs_id='rs1359790', effect_allele='A', beta_value=-0.077),
                        beta(rs_id='rs1470579', effect_allele='C', beta_value=0.103),
                        beta(rs_id='rs1552224', effect_allele='C', beta_value=-0.101),
                        beta(rs_id='rs17168486', effect_allele='T', beta_value=0.093),
                        beta(rs_id='rs1801282', effect_allele='G', beta_value=-0.113),
                        beta(rs_id='rs2237892', effect_allele='T', beta_value=-0.139),
                        beta(rs_id='rs2421016', effect_allele='T', beta_value=-0.056),
                        beta(rs_id='rs2796441', effect_allele='A', beta_value=-0.064),
                        beta(rs_id='rs2943640', effect_allele='C', beta_value=0.082),
                        beta(rs_id='rs3130501', effect_allele='G', beta_value=0.069),
                        beta(rs_id='rs340874', effect_allele='C', beta_value=0.060),
                        beta(rs_id='rs3794991', effect_allele='T', beta_value=0.130),
                        beta(rs_id='rs4457053', effect_allele='A', beta_value=-0.080),
                        beta(rs_id='rs459193', effect_allele='G', beta_value=0.078),
                        beta(rs_id='rs4607103', effect_allele='T', beta_value=-0.065),
                        beta(rs_id='rs4689388', effect_allele='A', beta_value=0.084),
                        beta(rs_id='rs516946', effect_allele='C', beta_value=0.078),
                        beta(rs_id='rs5215', effect_allele='T', beta_value=-0.075),
                        beta(rs_id='rs6813195', effect_allele='T', beta_value=-0.064),
                        beta(rs_id='rs7178572', effect_allele='G', beta_value=0.062),
                        beta(rs_id='rs7202877', effect_allele='G', beta_value=-0.106),
                        beta(rs_id='rs7578597', effect_allele='C', beta_value=-0.138),
                        beta(rs_id='rs7607980', effect_allele='C', beta_value=-0.093),
                        beta(rs_id='rs7674212', effect_allele='T', beta_value=-0.063),
                        beta(rs_id='rs7754840', effect_allele='C', beta_value=0.123),
                        beta(rs_id='rs7903146', effect_allele='T', beta_value=0.329),
                        beta(rs_id='rs7957197', effect_allele='A', beta_value=-0.077),
                        beta(rs_id='rs8042680', effect_allele='A', beta_value=0.063),
                        beta(rs_id='rs8050136', effect_allele='A', beta_value=0.107),
                        beta(rs_id='rs864745', effect_allele='C', beta_value=-0.085),
                        beta(rs_id='rs9505118', effect_allele='G', beta_value=-0.063)
                    ],
                    function_to_apply='add'
                )
            }
    }
)
if __name__ == '__main__':
    import vcf
    import sys
    from polygenic import polygenic
    allele_frequency_in_population = {
        'rs9847240': {'A': 0.2, 'C': 0.3, 'G': 0.4, 'T': 0.1},
        'rs2025905': {'A': 0.25, 'C': 0.35, 'G': 0.3, 'T': 0.1},
        'rs1835873': {'A': 0.15, 'C': 0.3, 'G': 0.45, 'T': 0.1},
        'rs1805005': {'A': 0.8, 'C': 0.05, 'G': 0.05, 'T': 0.1}
    }  # TODO make sure this dictionary covers all snips and sum == 1
    # model = HairColourModel(allele_frequency_in_population)
    data = {'rs10811661':'C/C',
         'rs10830963':'A/A',
         'rs10842994':'A/A',
         'rs11063069':'A/A',
         'rs1111875':'A/A',
         'rs11257655':'A/A',
         'rs11708067':'A/A',
         'rs12571751':'A/A',
         'rs12970134':'A/A',
         'rs13266634':'A/A',
         'rs1359790':'A/A',
         'rs1470579':'A/A',
         'rs1552224':'A/A',
         'rs17168486':'A/A',
         'rs1801282':'A/A',
         'rs2237892':'A/A',
         'rs2421016':'A/A',
         'rs2796441':'A/A',
         'rs2943640':'A/A',
         'rs3130501':'A/A',
         'rs340874':'A/A',
         'rs3794991':'A/A',
         'rs4457053':'A/A',
         'rs459193':'A/A',
         'rs4607103':'A/A',
         'rs4689388':'A/A',
         'rs516946':'A/A',
         'rs5215':'A/A',
         'rs6813195':'A/A',
         'rs7178572':'A/A',
         'rs7202877':'A/A',
         'rs7578597':'A/A',
         'rs7607980':'A/A',
         'rs7674212':'A/A',
         'rs7754840':'A/A',
         'rs7903146':'A/A',
         'rs7957197':'A/A',
         'rs8042680':'A/A',
         'rs8050136':'A/A',
         'rs864745':'A/A',
         'rs9505118':'A/A'}
    print(type_2_diabetes_model(data, allele_frequency_in_population))