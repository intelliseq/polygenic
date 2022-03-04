#!/bin/bash

# pharmvar.sh - a script to download pharmvar database

# url of the pharmvar database
blob:https://www.pharmvar.org/f7190663-b4cb-46ce-b943-134e91132877
url="ftp://ftp.ebi.ac.uk/pub/databases/pharmgkb/pharmvar/latest/pharmvar.tsv.gz"
url="ftp://ftp.ebi.ac.uk/pub/databases/pharmgkb/pharmgkb_latest.tar.gz"



# array of genes to download
export GENES=(


"CYP2A13", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2W1", "CYP3A4", "CYP3A5", "CYP3A7", "CYP3A43", "CYP4F2", "SLCO1B1", "DPYD", "NUDT15")

echo "haplotype_model:" > cyp2d6-pharmvar.yml
echo "  variants:" >> cyp2d6-pharmvar.yml
echo "    CYP2D6*1.001:" >> cyp2d6-pharmvar.yml
ls *.vcf | \
  sort -k1,1V | \
  xargs -i bash -c \
    'HAPLOTYPE=$(echo {} | sed "s|.vcf||" | sed "s|_|\*|"); echo "    $HAPLOTYPE:"; cat {} | grep -v "#" | sed "s|chr||"' | \
  awk '{if ($1 ~ /^22/) {chr=$1; pos=$2; ref=$4; alt=$5; id=chr"-"pos"-"ref"-"alt; print "      "id": {ref: \""ref"\", alt: \""alt"\", effect_allele: \""alt"\"}"} else {print $0}}' >> cyp2d6-pharmvar.yml