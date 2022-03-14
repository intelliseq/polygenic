#!/bin/bash
export VERSION=$(ls pharmvar-* | head -1 | grep -oP '[0-9]+\.[0-9]+\.[0-9]+')
mkdir -p models

function generate_model {
  GENE=$1
  GENE_LOWERCASE=$(echo "$GENE" | tr '[:upper:]' '[:lower:]')
  echo $GENE
  echo $GENE_LOWERCASE
  echo "haplotype_model:" > models/$GENE_LOWERCASE-pharmvar.yml
  echo "  variants:" >> models/$GENE_LOWERCASE-pharmvar.yml
  echo "    $GENE*1.001:" >> models/$GENE_LOWERCASE-pharmvar.yml
  ls pharmvar-$VERSION/$GENE/GRCh38/*.vcf | \
    sort -k1,1V | \
    xargs -i bash -c \
      'HAPLOTYPE=$(echo {} | sed "s|'pharmvar-$VERSION/$GENE/GRCh38/'||" | sed "s|.vcf||" | sed "s|_|\*|"); echo "    $HAPLOTYPE:"; cat {} | grep -v "#" | sed "s|chr||"' | \
  awk '{if ($1 ~ /^22/) {chr=$1; pos=$2; ref=$4; alt=$5; id=chr"-"pos"-"ref"-"alt; print "      "id": {ref: \""ref"\", alt: \""alt"\", effect_allele: \""alt"\"}"} else {print $0}}' >> models/$GENE_LOWERCASE-pharmvar.yml
}
export -f generate_model

ls pharmvar-$VERSION | xargs -i bash -c 'generate_model {}'
