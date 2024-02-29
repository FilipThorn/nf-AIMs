#!/bin/bash -l 

module load bioinfo-tools Nextflow/22.10.1 python/3.10.8 samtools/1.17 gcc/9.3.0 cairo/1.17.4 texinfo/6.8 texlive/2021 libcurl/7.45.0 readline/6.2-11 libicu/5.2-4 R/4.1.1 R_packages/4.1.1 java/OracleJDK_11.0.9 vcftools/0.1.16 biopython/1.80-py3.10.8 

NXF_HOME=/crex/proj/snic2020-16-126/Bop/MS2/08.AIMs/nf-AncestryPainter/

export NXF_DEFAULT_DSL=1

nextflow run main.nf -config profile/rackham.config -resume
