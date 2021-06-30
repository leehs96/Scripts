#!/bin/bash

# STAR \
#     --runThreadN 5 \
#     --runMode genomeGenerate \
#     --genomeDir /users/data/reference/ens \
#     --sjdbGTFfile /users/data/reference/ens/Homo_sapiens.GRCh38.102.chr.gtf \
#     --genomeFastaFiles users/data/reference/ens/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#     --sjdbOverhang 100

#rsem-prepare-reference --gtf /users/data/reference/ens/Homo_sapiens.GRCh38.102.chr.gtf \
#                        --star \
#                        --star-path ~/anaconda3/pkgs/star-2.7.6a-0/bin \
#                        -p 8 \
#                        /users/data/reference/ens/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                        /users/data/reference/ens/Homo_sapiens.GRCh38.dna.primary_assembly

coldata=$1

while IFS= read -r line || [ -n "$line" ] 
do
Line=$line

echo "$line"
FastqPath=`echo ${Line} | awk '{print $2}'`
Sample=`echo ${Line} | awk '{print $3}'`
OutputPath=`echo ${Line} | awk '{print $4}'`
 
if [ ! -d ${OutputPath} ]; then
  mkdir ${OutputPath}
fi

mkdir ${OutputPath}/trim
mkdir ${OutputPath}/bam
mkdir ${OutputPath}/RSEM
mkdir ${OutputPath}/featureCounts


echo ${Line}

trimmomatic PE \
    -threads 5 \
    ${FastqPath}/${Sample}_RNA_1.fastq.gz \
    ${FastqPath}/${Sample}_RNA_2.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_U_1.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_U_2.fastq.gz \
    ILLUMINACLIP:/users/data/reference/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:5 \
    TRAILING:5 \
    SLIDINGWINDOW:5:10 \
    MINLEN:20 && 

rm  ${OutputPath}/trim/${Sample}_RNA_U_1.fastq.gz 
rm  ${OutputPath}/trim/${Sample}_RNA_U_2.fastq.gz 

wait

STAR \
   --genomeDir \
   /users/data/reference/ens \
   --runThreadN 6 \
   --sjdbGTFfile /users/data/reference/ens/Homo_sapiens.GRCh38.102.chr.gtf \
   --readFilesIn ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz \
   ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz \
   --sjdbOverhang 100 \
   --readFilesCommand zcat \
   --outSAMtype BAM SortedByCoordinate \
   --outFilterMultimapNmax 20 \
   --outFilterType BySJout \
   --twopassMode Basic \
   --quantMode TranscriptomeSAM GeneCounts \
   --outFileNamePrefix ${OutputPath}/bam/${Sample}.STAR. &&

wait 

samtools index ${OutputPath}/bam/${Sample}.STAR.Aligned.sortedByCoord.out.bam &&

UnS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $2}'`
FS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $3}'`
RS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $4}'`

diff = $((FS - RS))

if [ ${diff#-} -gt $(((FS + RS) / 2)) ]
then
    if [ ${FS} > ${RS} ]
    then
        dicision='Stranded'
    elif [ ${FS} < ${RS} ]
    then
        dicision='ReverselyStranded'
    elif [ ${FS} == ${RS} ]
    then
        dicision='UnStranded'
else
    dicision='UnStranded'
fi


if [ ${dicision} == 'UnStranded' ]
then

    RSEM_Strandedness_dicision= 0.5
    featureCounts_Strandedness_dicision= 0

elif [ ${dicision} == 'Stranded' ]
then

    RSEM_Strandedness_dicision= 1
    featureCounts_Strandedness_dicision= 1

elif [ ${dicision} == 'ReverselyStranded' ]
then

    RSEM_Strandedness_dicision= 0
    featureCounts_Strandedness_dicision= 2

else
    echo "ERROR undefined RSEM strandedness"
    exit 1


rsem-calculate-expression \
    -p 6 \
    --alignments \
    --paired-end \
    --forward-prob ${RSEM_Strandedness_dicision} \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/ens/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 

wait

rm ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam
rm  ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz
rm  ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz 

done < ${coldata}

wait

featureCounts -T 10 -s ${featureCounts_Strandedness_dicision} -p \
    -a /users/data/reference/ens/Homo_sapiens.GRCh38.102.chr.gtf \
    -o ${OutputPath}/featureCounts/featurecounts.results.txt \
    ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam &&

sed -i "s:${OutputPath}/bam/::" ${OutputPath}/featureCounts/featurecounts.results.txt
sed -i "s:.STAR.Aligned.sortedByCoord.out.bam::" ${OutputPath}/featureCounts/featurecounts.results.txt

wait

cut -f1,7- ${OutputPath}/featureCounts/featurecounts.results.txt > ${OutputPath}/featureCounts/featurecounts.results.final.txt

echo "All DONE"