#!/bin/bash
rm test.sample.vcf.*
cat test.sample.vcf | bgzip -c > test.sample.vcf.gz
tabix -p vcf test.sample.vcf.gz
