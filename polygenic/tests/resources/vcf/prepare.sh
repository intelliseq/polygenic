#!/bin/bash
rm mini.sample.vcf.*
cat mini.sample.vcf | bgzip -c > mini.sample.vcf.gz
tabix -p vcf mini.sample.vcf.gz
