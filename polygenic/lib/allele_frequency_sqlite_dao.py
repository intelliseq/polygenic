# import sqlite3
# import os
# import sys

# module_path = os.path.abspath(__file__).rsplit(os.path.sep, 5)[0]
# sys.path.insert(0, module_path)

# import argparse
# from collections import OrderedDict
# from typing import List
# from polygenic.lib.gnomad_parser import SNPdata
# from polygenic.lib.gnomad_parser import get_fields
# from polygenic.lib.gnomad_parser import get_data_with_unique_rsid_and_non_unique_containing_most_frequent_alt

# DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)).rsplit('/', 2)[0], 'resources', 'databases',
#                        'allele_frequency_db_gnomad_2.1.1.sqlite')#'allele_frequency_db_gnomad_3.0.sqlite')


# def create_table(path: str = DB_PATH):
#     if not os.path.exists(os.path.dirname(DB_PATH)):
#         os.makedirs(os.path.dirname(DB_PATH))
#     with sqlite3.connect(path) as conn:
#         cur = conn.cursor()
#         # query = '''
#         #             CREATE TABLE IF NOT EXISTS allele_frequency
#         #                 (
#         #                     rsid VARCHAR PRIMARY KEY,
#         #                     ref VARCHAR,
#         #                     alt VARCHAR,
#         #                     AF FLOAT,
#         #                     AF_nfe FLOAT,
#         #                     AF_eas FLOAT,
#         #                     AF_afr FLOAT,
#         #                     AF_amr FLOAT,
#         #                     AF_asj FLOAT,
#         #                     AF_fin FLOAT,
#         #                     AF_ami FLOAT,
#         #                     AF_oth FLOAT,
#         #                     id_unique BOOLEAN
#         #                 )
#         #         '''
#         query = '''
#                     CREATE TABLE IF NOT EXISTS allele_frequency
#                         (
#                             rsid VARCHAR PRIMARY KEY,
#                             ref VARCHAR,
#                             alt VARCHAR,
#                             AF FLOAT,
#                             AF_nfe FLOAT,
#                             AF_eas FLOAT,
#                             AF_afr FLOAT,
#                             AF_amr FLOAT,
#                             AF_asj FLOAT,
#                             AF_fin FLOAT,
#                             AF_oth FLOAT,
#                             id_unique BOOLEAN
#                         )
#                 '''
#         cur.execute(query)


# def insert_many(data: List[SNPdata], path: str = DB_PATH):
#     with sqlite3.connect(path) as conn:
#         cur = conn.cursor()
#         #query = 'INSERT INTO allele_frequency VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'
#         query = 'INSERT INTO allele_frequency VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'
#         cur.executemany(query, data)


# def read_multiple_rs(rs_ids: List[str], path: str = DB_PATH) -> List[SNPdata]:
#     if not os.path.exists(path):
#         raise IOError('Path {} does not exist'.format(path))
#     with sqlite3.connect(path) as conn:
#         cur = conn.cursor()
#         query = 'SELECT * FROM allele_frequency WHERE rsid IN ({})'.format(
#             ', '.join('"{}"'.format(rs_id) for rs_id in rs_ids))
#         cur.execute(query)
#         return [SNPdata(*x) for x in cur.fetchall()]


# def get_column_names(path: str = DB_PATH) -> List[str]:
#     with sqlite3.connect(path) as conn:
#         cur = conn.execute('select * from allele_frequency')
#         return [description[0] for description in cur.description]


# if __name__ == '__main__':
#     #print(get_column_names())
#     parser = argparse.ArgumentParser(description='')  # todo dodać opis
#     parser.add_argument('data_dir', type=str, help='Directory containing gnomad data')
#     parsed_args = parser.parse_args()
#     create_table()
#     data = []
#     # for snpdata in get_fields('/home/wojtek/Pobrane/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz', ['ID', 'REF', 'ALT'],
#     #                           OrderedDict([('AF', float), ('AF_nfe', float), ('AF_eas', float), ('AF_afr', float),
#     #                                        ('AF_amr', float), ('AF_asj', float), ('AF_fin', float), ('AF_oth', float),
#     #                                        ('AF_sas', float)])):
#     #     if snpdata.rsid:
#     #         data.append(snpdata)
#     # print('filtering')
#     # data = get_data_with_unique_rsid_and_non_unique_containing_most_frequent_alt(data)
#     # print(len(data))
#     # print(data[-1])
#     # insert_many(data)
#     #raise

#     for chrom in list(range(1, 23)) + ['X', 'Y']:
#         print("Processing chromosome {}".format(chrom))
#         file_name = 'gnomad.genomes.r3.0.sites.chr{}.vcf.bgz'.format(chrom)
#         file_path = os.path.join(parsed_args.data_dir, file_name)
#         data = []
#         for snpdata in get_fields(file_path, ['ID', 'REF', 'ALT'],
#                                   OrderedDict([('AF', float), ('AF_nfe', float), ('AF_eas', float), ('AF_afr', float),
#                                                ('AF_amr', float), ('AF_asj', float), ('AF_fin', float),
#                                                ('AF_ami', float), ('AF_oth', float)])):
#             if snpdata.rsid and snpdata.rsid != '.':
#                 data.append(snpdata)
#         print('filtering')
#         data = get_data_with_unique_rsid_and_non_unique_containing_most_frequent_alt(data)
#         print("inserting data")
#         insert_many(data)
#     print(read_multiple_rs(['rs59514816', 'rs1214589876']))
