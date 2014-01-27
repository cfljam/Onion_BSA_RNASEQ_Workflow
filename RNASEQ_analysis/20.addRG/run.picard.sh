cat file.lst | while read filename outname; 
do 
    java -jar /software/x86_64/picard-tools-1.79/AddOrReplaceReadGroups.jar input=$filename output=${outname}_RG.bam RGID=${outname} RGLB=PE RGPL=illumina RGPU=barcode RGSM=${outname}; 
done