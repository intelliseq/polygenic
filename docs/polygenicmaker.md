# polygenicmaker

## Preparing models from Global Biobank Engine
```
OUTPUT_DIRECTORY=/tmp/output
docker run -v "$OUTPUT_DIRECTORY":/output marpiech/polygenicmaker:2.0.25 polygenicmaker gbe-model --code HC710
```

## Preparing models from Pan Biobank UK
#### Downloading index
```
polygenicmaker \
  biobankuk-index \
  --output-directory /tmp/polygenic/results
```
#### Preparing model
```
polygenicmaker \
  biobankuk-model \
  --code 2395  \
  --sex both_sexes \
  --coding 4 \
  --index /tmp/polygenic/results/panukbb_phenotype_manifest.tsv \
  --output-directory /tmp/polygenic/results/model \
  --variant-metrics /tmp/polygenic/results/full_variant_qc_metrics.txt \
  --threshold 1e-08 \
  --source-ref-vcf /tmp/polygenic/results/ALL.2of4intersection.20100804.genotypes.vcf.gz \
  --target-ref-vcf /tmp/marpiech/kenobi/resources/GRCh38.dbSNP155.chr.norm.rsidonly.vcf.gz
```
