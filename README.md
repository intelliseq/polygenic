# polygenic

[![PyPI](https://img.shields.io/pypi/v/polygenic.svg)](https://pypi.python.org/pypi/polygenic)

python package for computation of polygenic scores based for particular sample

## How to install
### Using pip
```
pip3 install --upgrade polygenic
```
### In new conda environment
```
docker run -it conda/miniconda3 /bin/bash
```
```
yes | conda create --name py38 python=3.8
eval "$(conda shell.bash hook)"
conda activate py38
### should be 3.8
python --version

### gcc is missing to build pytabix
apt -qq update
apt -y install build-essential

pip install polygenic
```

## How to run
```
polygenic --vcf [your_vcf_gz] --model [your_model] [other raguments]
```

### Arguments
#### Required
- `--vcf` vcf.gz file with genotypes (tabix index should be available)
- `--model` path to model file
#### Optional
- `--log_file` log file
- `--out_dir` directory for result jsons
- `--population` population code
- `--models_path` path to a directory containing models
- `--af` an indexed vcf.gz file containing allele freq data
- `--version` prints version of package

## Building models
Models are pure python scripts tha use "sequencing query languange" called seqql.  
It is required to import language elements.
```
from polygenic.lib.model.seqql import PolygenicRiskScore
from polygenic.lib.model.seqql import ModelData
from polygenic.lib.model.category import QuantitativeCategory
```

It recommended to add variable pointing population for which score was prepared
```
trait_was_prepared_for_population = "eas"
```
The list of accepted population identifiers:
- `nfe` - Non-Finnish European ancestry
- `eas` - East Asian ancestry
- `afr` - African-American/African ancestry
- `amr` - Latino ancestry
- `asj` - Ashkenazi Jewish ancestry,
- `fin` - Finnish ancestry
- `oth` - Other ancestry

The most important part of model is model itself. Currently it is possible to use PolygenicRiskScore
```
model = PolygenicRiskScore(categories = ..., snps_and_coeffcients = ..., model_type = ...)
```

categories is a list of named results ranges (`QuantitativeCategory`) that can be used to define bucket for which interpretation will be generated
```QuantitativeCategory(from_= ..., to=..., category_name=...)```

snps_and_coeffcients is a list of snps
with their effect allele in genomic notation and coeffcient value. Snps are defined by their rsid
```
'rs10012': ModelData(effect_allele='G', coeff_value=0.369215857410143),
```
## Example model
```
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
```

## Rescaling model results
It is possible to further rescale model results within each Category
```
categories=[
        QuantitativeCategory(from_=1.371624087, to=2.581880425, category_name='High risk', scale_from = 2, scale_to = 3),
        QuantitativeCategory(from_=1.169616034, to=1.371624087, category_name='Potential risk', scale_from = 1, scale_to = 2),
        QuantitativeCategory(from_=-0.346748358, to=1.169616034, category_name='Average risk', scale_from = 0, scale_to = 1),
	    QuantitativeCategory(from_=-1.657132197, to=-0.346748358, category_name='Low risk', scale_from = -1, scale_to = 0)
    ],
```
