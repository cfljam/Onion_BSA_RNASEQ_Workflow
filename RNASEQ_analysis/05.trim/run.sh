#!/bin/sh
echo "trim fastq files"
cat sample.lst | while read line; do fq2trimmed_extended.pl -c=trim -5=15 -3=5 -f1=D1B7JACXX_NZGL00123_${line}_R1_001.fastq -f2=D1B7JACXX_NZGL00123_${line}_R2_001.fastq -o1=${line}_q15_p15_R1.fastq -o2=${line}_q15_p15_R2.fastq; done
