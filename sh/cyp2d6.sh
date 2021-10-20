#!/bin/bash
echo "haplotype_model:" > cyp2d6-pharmvar.yml
echo "  variants:" >> cyp2d6-pharmvar.yml
echo "    CYP2D6*1.001:" >> cyp2d6-pharmvar.yml
ls *.vcf | \
  sort -k1,1V | \
  xargs -i bash -c \
    'HAPLOTYPE=$(echo {} | sed "s|.vcf||" | sed "s|_|\*|"); echo "    $HAPLOTYPE:"; cat {} | grep -v "#" | sed "s|chr||"' | \
  awk '{if ($1 ~ /^22/) {chr=$1; pos=$2; ref=$4; alt=$5; id=chr"-"pos"-"ref"-"alt; print "      "id": {ref: \""ref"\", alt: \""alt"\", effect_allele: \""alt"\"}"} else {print $0}}' >> cyp2d6-pharmvar.yml