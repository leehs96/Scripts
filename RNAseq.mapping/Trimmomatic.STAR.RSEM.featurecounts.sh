#!/bin/bash

###################################################################
#                                                                 #
#  STAR version : 2.7.4a                                          #
#                                                                 #
#  RSEM version : 1.3.1                                           #
#                                                                 #
#  reference fasta : Homo_sapiens.GRCh38.dna.primary_assembly.fa  #
#                                                                 #
#  reference GTF : Homo_sapiens.GRCh38.104.gtf                    #
#                                                                 #    
###################################################################

# in case of reference file malfunction

#STAR \
#     --runThreadN 5 \
#     --runMode genomeGenerate \
#     --genomeDir /users/data/reference/hg38/ \
#     --sjdbGTFfile /users/data/reference/hg38/Homo_sapiens.GRCh38.104.chr.gtf \
#     --genomeFastaFiles users/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#     --sjdbOverhang 100

# rsem-prepare-reference --gtf /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
#                        --star \
#                        --star-path ~/anaconda3/pkgs/star-2.7.6a-0/bin \
#                        -p 8 \
#                        /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                        /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly

coldata=$1
dofeaturecount=$2  #y = do featurecounts, n = don't

while IFS= read -r line || [ -n "$line" ] 
do
Line=$line

echo "$line"
FastqPath=`echo ${Line} | awk '{print $2}'`
Sample=`echo ${Line} | awk '{print $3}'`
OutputPath=`echo ${Line} | awk '{print $4}'`
 
if [ ! -d ${OutputPath} ]
then
  mkdir ${OutputPath}
fi

dir=("trim" "bam" "RSEM" "FeatureCounts")

for dirname in ${dir[@]}
do
if [ ! -d ${OutputPath}/${dirname} ]
then
  mkdir ${OutputPath}/${dirname}
fi
done

rm /users/data/log/${Sample}.log

date +"%d-%m-%Y %T: ${Sample} RNA processing START" |tee -a /users/data/log/${Sample}.log

echo -e "${Sample} : Outputpath = ${OutputPath}  " |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : trimommatic  " |tee -a /users/data/log/${Sample}.log

trimmomatic PE \
    -threads 5 \
    ${FastqPath}/${Sample}_RNA_1.fastq.gz \
    ${FastqPath}/${Sample}_RNA_2.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_U_1.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz \
    ${OutputPath}/trim/${Sample}_RNA_U_2.fastq.gz \
    ILLUMINACLIP:/users/data/reference/hg38/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:5 \
    TRAILING:5 \
    SLIDINGWINDOW:5:10 \
    MINLEN:20 && 

rm  ${OutputPath}/trim/${Sample}_RNA_U_1.fastq.gz 
rm  ${OutputPath}/trim/${Sample}_RNA_U_2.fastq.gz 

echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : trimommatic Done \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log

wait
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : STAR \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
STAR \
   --genomeDir \
   /users/data/reference/hg38/ \
   --runThreadN 6 \
   --sjdbGTFfile /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
   --readFilesIn ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz \
   ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz \
   --readFilesCommand zcat \
   --outSAMtype BAM SortedByCoordinate \
   --outFilterMultimapNmax 20 \
   --outFilterType BySJout \
   --twopassMode Basic \
   --quantMode TranscriptomeSAM GeneCounts \
   --outFileNamePrefix ${OutputPath}/bam/${Sample}.STAR. &&

wait 

echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : STAR Done \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
samtools index ${OutputPath}/bam/${Sample}.STAR.Aligned.sortedByCoord.out.bam &&

UnS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $2}'`
FS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $3}'`
RS=`grep 'ENSG00000111640' ${OutputPath}/bam/${Sample}.STAR.ReadsPerGene.out.tab | awk '{print $4}'`

distance=$((FS - RS))

if [ ${distance#-} -gt $(((FS + RS) / 2)) ]
then
    if [ ${FS} -gt ${RS} ]
    then
        dicision='Stranded'
    elif [ ${FS} -lt ${RS} ]
    then
        dicision='ReverselyStranded'
    elif [ ${FS} -eq ${RS} ]
    then
        dicision='UnStranded'
    else
        echo "ERROR undefined type of strandedness"
        exit 1
    fi
else
    dicision='UnStranded'
fi



echo -e "------------------------------\n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : ${dicision}" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------\n\n" |tee -a /users/data/log/${Sample}.log


echo -e "------------------------------\n\n" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : RSEM \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------\n\n" |tee -a /users/data/log/${Sample}.log

if [ ${dicision} = "UnStranded" ]
then

    rsem-calculate-expression \
    -p 6 \
    --alignments \
    --paired-end \
    --forward-prob 0.5 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample}


elif [ ${dicision} = "Stranded" ]
then

    rsem-calculate-expression \
    -p 6 \
    --alignments \
    --paired-end \
    --forward-prob 1 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 

elif [ ${dicision} = "ReverselyStranded" ]
then

    rsem-calculate-expression \
    -p 6 \
    --alignments \
    --paired-end \
    --forward-prob 0 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 

else
    echo "ERROR undefined RSEM strandedness"
    exit 1
fi &&

rm ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam
rm ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz
rm ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz 

echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : RSEM Done \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log

wait


done < ${coldata}


wait
if [ ${dofeaturecount} = "y" ]
then

echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : FeatureCounts \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log


    if [ ${dicision} = "UnStranded" ]
    then

        featureCounts -T 10 -s 0 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 



    elif [ ${dicision} = "Stranded" ]
    then

        featureCounts -T 10 -s 1 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 


    elif [ ${dicision} = "ReverselyStranded" ]
    then

        featureCounts -T 10 -s 2 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 


    else
        echo "ERROR undefined RSEM strandedness"
        exit 1
    fi &&

    sed -i "s:${OutputPath}/bam/::" ${OutputPath}/FeatureCounts/featurecounts.results.txt
    sed -i "s:.STAR.Aligned.sortedByCoord.out.bam::" ${OutputPath}/featureCounts/featurecounts.results.txt
    cut -f1,7- ${OutputPath}/FeatureCounts/featurecounts.results.txt > ${OutputPath}/featureCounts/featurecounts.results.final.txt
elif [ ! ${dofeaturecount} = "n" ]
then
    echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
    date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
    echo -e "${Sample} : All Done w/o FeatureCounts \n\n" |tee -a /users/data/log/${Sample}.log
    echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
else
    echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
    date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
    echo -e "${Sample} : 1 = do featurecounts, 0 = don't \n\n" |tee -a /users/data/log/${Sample}.log
    echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
    exit
fi &&

echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample}" |tee -a /users/data/log/${Sample}.log
echo -e "${Sample} : All Done w/ FeatureCounts \n\n" |tee -a /users/data/log/${Sample}.log
echo -e "------------------------------ \n\n" |tee -a /users/data/log/${Sample}.log

#fin.