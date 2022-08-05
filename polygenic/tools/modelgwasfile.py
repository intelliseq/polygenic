from polygenic.data.gwas_data import GwasData as Gwas

def run(args):
    gwas = Gwas()
    gwas.load_gwas_data_from_csv(args.gwas_file, vars(args))

    # dbsnp37 = vcf_accessor.VcfAccessor("/home/marpiech/data/vcf/dbsnp155.grch37.norm.vcf.gz")
    # dbsnp38 = vcf_accessor.VcfAccessor("/home/marpiech/data/vcf/dbsnp155.grch38.norm.vcf.gz")

    # data = csv.get_data()

    # print(str(data.iloc[0:2, :]))