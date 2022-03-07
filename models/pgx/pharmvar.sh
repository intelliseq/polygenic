#!/bin/bash

# pharmvar.sh - a script to download pharmvar database

# script directory
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# working directory
export WORKING_DIR=/tmp

# url of the pharmvar database
# curl 'https://www.pharmvar.org/get-download-file?name=ALL&refSeq=ALL&fileType=zip&version=current' \
#   -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:97.0) Gecko/20100101 Firefox/97.0' \
#   -H 'Accept: */*' \
#   -H 'Accept-Language: en-US,en;q=0.5' \
#   -H 'Accept-Encoding: gzip, deflate, br' \
#   -H 'Connection: keep-alive' \
#   -H 'Referer: https://www.pharmvar.org/download' \
#   -H 'Sec-Fetch-Dest: empty' \
#   -H 'Sec-Fetch-Mode: cors' \
#   -H 'Sec-Fetch-Site: same-origin' \
#   -H 'DNT: 1' \
#   -H 'Sec-GPC: 1' \
#   --output $WORKING_DIR/pharmvar.zip

# unzip $WORKING_DIR/pharmvar.zip -d /tmp/

export VERSION=$(ls $WORKING_DIR | grep "pharmvar-" | sed 's/.*-//' | sort -r | head -n 1)

function parse_variant {

    local VCF=$1
    local GENE=$2
    local HAPLOTYPE=$GENE"*"$(echo $VCF | sed 's/.vcf//' | sed 's/.*_//' | sed 's/.*rs/rs/')
    echo "    $HAPLOTYPE:"; 
    cat $VCF | grep -v '#' | sed "s|chr||" | \
    awk '{chr=$1; pos=$2; ref=$4; alt=$5; id=chr"-"pos"-"ref"-"alt; print "      "id": {ref: \""ref"\", alt: \""alt"\", effect_allele: \""alt"\"}"}'

}
export -f parse_variant

function parse_gene {
    local GENE=$1
    local GENE_LOWERCASE=$(echo $GENE | tr '[:upper:]' '[:lower:]')
    local MODEL_FILE="$GENE_LOWERCASE-pharmvar-$VERSION.yml"
    echo "haplotype_model:" > $MODEL_FILE
    echo "  variants:" >> $MODEL_FILE
    echo "    $GENE*1.001:" >> $MODEL_FILE
    ls /tmp/pharmvar-$VERSION/$GENE/GRCh38/*.vcf | \
      sort -k1,1V | \
      xargs -i bash -c "parse_variant {} $GENE >> $MODEL_FILE"
}
export -f parse_gene

ls /tmp/pharmvar-$VERSION/ | xargs -i bash -c "parse_gene {}"


