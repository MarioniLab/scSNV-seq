samtools view -hb -r $arg PATH_TO_TAPESTRI_OUTPUT/results/bam/TAPESTRI_OUTPUT_NAME.tube1.cells.bam > tmp/$arg.bam

samtools sort -o tmp/$arg.sorted.bam tmp/$arg.bam

samtools index tmp/$arg.sorted.bam

gatk HaplotypeCaller --java-options "-Xmx4g" -R PATH_TO_REFERENCE_WITH_pKLV2/hg38_pKLV2_bwa.fa  -I  tmp/$arg.sorted.bam -O vcf/$arg.vcf.gz -ploidy 3

rm tmp/$arg.sorted.bam
rm tmp/$arg.bam
rm rm tmp/$arg.sorted.bam.bai

# usage per cell: arg is the cell barcode
# hg38_pKLV2_bwa.fa is the hg38 reference genome with added pKLV2 (for the genotype calling only (without barcodes) adding pKLV2 is not needed)


