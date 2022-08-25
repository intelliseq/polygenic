#!/bin/bash
rm test-vcf-general.vcf.*
cat test-vcf-general.vcf | bgzip -c > test-vcf-general.vcf.gz
tabix -p vcf test-vcf-general.vcf.gz
