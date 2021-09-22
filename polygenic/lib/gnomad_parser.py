# import gzip
# import re
# import numpy as np
# from collections import OrderedDict
# from typing import List
# from typing import Dict
# from typing import Callable
# from collections import namedtuple

# ##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">
# ##INFO=<ID=AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Non-Finnish European ancestry">
# ##INFO=<ID=AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry">
# ##INFO=<ID=AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American/African ancestry">
# ##INFO=<ID=AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry">
# ##INFO=<ID=AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry">
# ##INFO=<ID=AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry">
# ##INFO=<ID=AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of Other ancestry">

# #SNPdata = namedlist('SNPdata', 'rsid ref alt AF AF_nfe AF_eas AF_afr AF_amr AF_asj AF_fin AF_ami AF_oth is_unique',
# SNPdata = namedtuple('SNPdata', ['rsid', 'ref', 'alt', 'AF', 'AF_nfe', 'AF_eas', 'AF_afr', 'AF_amr', 'AF_asj', 'AF_fin', 'AF_oth', 'is_unique'])


# def get_fields(path: str, basic_fields: List[str], info_fields: Dict[str, Callable]):
#     with gzip.open(path, 'rt') as f:
#         for line in f:
#             if not line.startswith('#'):
#                 CHROM, POS, ID, REF, ALT, QUAL, FILTER, info = line.split()
#                 d = {'CHROM':CHROM, 'POS':POS, 'ID':ID, 'REF':REF, 'ALT':ALT, 'QUAL':QUAL, 'FILTER':FILTER}
#                 out_dict = {x: d[x] for x in basic_fields}
#                 out_dict.update(get_fields_from_info_string(info, info_fields))
#                 fields = basic_fields + list(info_fields.keys())
#                 yield SNPdata(*([out_dict[field] for field in fields]))


# def get_fields_from_info_string(info: str, fields: Dict[str, Callable]):
#     ret_dict = {}
#     for field, type_ in fields.items():
#         found = re.search('{}=(.+?);'.format(field), info)
#         try:
#             ret_dict[field] = type_(found.group(1))
#         except AttributeError:
#             ret_dict[field] = None
#     return ret_dict


# def get_data_with_unique_rsid_and_non_unique_containing_most_frequent_alt(data: List[SNPdata]):
#     data.sort(key=lambda snp_data: snp_data.rsid)
#     unique_data = []
#     non_unique_data = {}
#     i = 0
#     while i < len(data):
#         j = i + 1
#         while j < len(data) and data[i].rsid == data[j].rsid:
#             append_data_to_list_value_of_dict(non_unique_data, data[i].rsid, data[j])
#             j += 1
#         if j == i + 1:
#             unique_data.append(data[i])
#         else:
#             append_data_to_list_value_of_dict(non_unique_data, data[i].rsid, data[i])
#         i = j
#     for element in unique_data:
#         element.is_unique = True
#     for v in non_unique_data.values():
#         data_containing_most_frequent_alt_allele = get_data_containing_most_frequent_alt_allele(v)
#         data_containing_most_frequent_alt_allele.is_unique = False
#         unique_data.append(data_containing_most_frequent_alt_allele)
#     return unique_data


# def append_data_to_list_value_of_dict(dictionary: Dict[str, List[SNPdata]], key: str, value_to_append: SNPdata):
#     if key not in dictionary:
#         dictionary[key] = [value_to_append]
#     else:
#         dictionary[key].append(value_to_append)


# def get_data_containing_most_frequent_alt_allele(data: List[SNPdata]):
#     return max(data, key=lambda snp_data: cast_nonetype_to_negative_inf(snp_data.AF))


# def cast_nonetype_to_negative_inf(float_or_nonetype):
#     if float_or_nonetype is None:
#         return np.NINF
#     return float_or_nonetype


# if __name__ == '__main__':
#     data = []
#     for snpdata in get_fields('/home/wojtek/Pobrane/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz', ['ID', 'REF', 'ALT'],
#                               OrderedDict([('AF', float), ('AF_nfe', float), ('AF_eas', float), ('AF_afr', float),
#                                           ('AF_amr', float), ('AF_asj', float), ('AF_fin', float), ('AF_oth', float)])):
#         if snpdata.rsid and snpdata.rsid != '.':
#             data.append(snpdata)
#     print('filtering')
#     data = get_data_with_unique_rsid_and_non_unique_containing_most_frequent_alt(data)
