import sys
import polygenic.tools as tools

from polygenic.error.polygenic_exception import PolygenicException

def main(args=sys.argv[1:]):



    try:
        if args[0] == 'pgs-compute':
            tools.pgscompute.main(args[1:])
        elif args[0] == 'model-biobankuk':
            tools.modelbiobankuk.main(args[1:])
        elif args[0] == 'model-pgscat':
            tools.modelpgscat.main(args[1:])
        elif args[0] == 'vcf-index':
            tools.vcfindex.main(args[1:])
        else:
            print('ERROR: Please select proper tool name"')
            print("""
            Program: polygenic toolkit (downloads gwas data, builds and computes polygenic scores)
            Contact: Marcin Piechota <piechota@intelliseq.com>
            Usage:   pgstk <command> [options]

            Commands:
            pgs-compute             computes pgs score for vcf file
            model-biobankuk         prepare polygenic score model based on gwas results from biobankuk
            model-pgscat            prepare polygenic score model based on gwas results from PGS Catalogue
            vcf-index               prepare rsidx for vcf
            plot-gwas               manhattan plot for gwas results
            """)
    except PolygenicException as e:
        print("Analysis failed")
        print("ERROR: " + str(e))
    except RuntimeError as e:
        print("ERROR: " + str(e))

if __name__ == '__main__':
    main(sys.argv[1:])
