cat lists.lst | while read line;
do samtools merge ${line}.bam ${line}_L002.bam ${line}_L003.bam;samtools sort ${line}.bam ${line}_s
done