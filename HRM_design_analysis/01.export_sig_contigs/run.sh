#!/bin/sh
echo "form fasta of significant contigs"
cut -f1 ../50.BoltVsNonCMH_V2/sig_hits | uniq | fastaselect -c -f ../../NZGL00123a/41.blast_mt_DNA/11.remove_mt_hits/filtered_assembly > BSA_hits.fasta
