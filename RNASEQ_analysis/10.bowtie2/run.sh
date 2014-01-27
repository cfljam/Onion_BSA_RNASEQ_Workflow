#!/bin/sh
echo "run bowtie2 on onion trimmed fastq files"
cat sample.lst | while read line; do bowtie2 --sensitive -X 500 --un ${line}_unaligned -S ${line}.sam -x {{db}} -1 ../05.trim/${line}_q15_p15_R1.fastq -2 ../05.trim/${line}_q15_p15_R2.fastq; done
