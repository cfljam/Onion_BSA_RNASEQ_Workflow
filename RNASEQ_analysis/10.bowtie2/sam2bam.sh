for x in *.sam;
do 
    name=$(basename "${x}" .sam); 
    samtools view -S -b ${x} > ${name}.bam; 
done