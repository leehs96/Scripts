#!/bin/bash
set -e
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

FastqPath=`echo ${Line} | awk '{print $2}'`
Sample=`echo ${Line} | awk '{print $3}'`
OutputPath=`echo ${Line} | awk '{print $4}'`

if [ ${dofeaturecount} = 'y' ]
then
echo -e "\n\n\nDo FeatureCounts after processing of all samples\n\n\n"
elif [ ${dofeaturecount} = 'n' ]
then
echo -e "\n\n\nDo NOT FeatureCounts after processing of all samples\n\n\n"
else
echo -e "\n\n\nERROR : please, assign 2nd parameter for featurecounts : y = do featurecounts, n = don't"
exit 1
fi

if [ ! -d ${OutputPath} ]
then
  mkdir ${OutputPath}
fi

dir=("trim" "bam" "RSEM" "FeatureCounts" "log")

for dirname in ${dir[@]}
do
if [ ! -d ${OutputPath}/${dirname} ]
then
  mkdir ${OutputPath}/${dirname}
fi
done

if [ -f ${OutputPath}/log/${Sample}.log ]
then
rm ${OutputPath}/log/${Sample}.log
fi

echo -e "---------------------------------------------------- \n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T: ${Sample} RNA processing START" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "\n---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log


echo -e "${Sample} : Fastq path = ${FastqPath}  " |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : Output path = ${OutputPath}  " |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : trimommatic  \n\n\n" |tee -a ${OutputPath}/log/${Sample}.log



trimmomatic PE \
    -threads 10 \
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
    MINLEN:20 2>> ${OutputPath}/log/${Sample}.log

rm  ${OutputPath}/trim/${Sample}_RNA_U_1.fastq.gz 
rm  ${OutputPath}/trim/${Sample}_RNA_U_2.fastq.gz 

echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : trimommatic Done \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log

wait
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : STAR \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
STAR \
   --genomeDir \
   /users/data/reference/hg38/ \
   --runThreadN 10 \
   --sjdbGTFfile /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
   --readFilesIn ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz \
   ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz \
   --readFilesCommand zcat \
   --outSAMtype BAM SortedByCoordinate \
   --outFilterMultimapNmax 20 \
   --outFilterType BySJout \
   --twopassMode Basic \
   --quantMode TranscriptomeSAM GeneCounts \
   --outFileNamePrefix ${OutputPath}/bam/${Sample}.STAR. 2>> ${OutputPath}/log/${Sample}.log

wait 

echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : STAR Done \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
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
        echo "ERROR undefined type of strandedness" |tee -a ${OutputPath}/log/${Sample}.log
        exit 1
    fi
else
    dicision='UnStranded'
fi



echo -e "----------------------------------------------------\n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample}'s RNA library sequenced as ${dicision}\n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "----------------------------------------------------\n\n" |tee -a ${OutputPath}/log/${Sample}.log


echo -e "----------------------------------------------------\n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : RSEM \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "----------------------------------------------------\n\n" |tee -a ${OutputPath}/log/${Sample}.log

if [ ${dicision} = "UnStranded" ]
then

    rsem-calculate-expression \
    -p 10 \
    --alignments \
    --paired-end \
    --forward-prob 0.5 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 2>> ${OutputPath}/log/${Sample}.log


elif [ ${dicision} = "Stranded" ]
then

    rsem-calculate-expression \
    -p 10 \
    --alignments \
    --paired-end \
    --forward-prob 1 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 2>> ${OutputPath}/log/${Sample}.log

elif [ ${dicision} = "ReverselyStranded" ]
then

    rsem-calculate-expression \
    -p 10 \
    --alignments \
    --paired-end \
    --forward-prob 0 \
    --bam --no-bam-output \
    ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam \
    /users/data/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly \
    ${OutputPath}/RSEM/${Sample} 2>> ${OutputPath}/log/${Sample}.log

else
    echo "ERROR : undefined RSEM strandedness" |tee -a ${OutputPath}/log/${Sample}.log
    exit 1
fi &&

rm ${OutputPath}/bam/${Sample}.STAR.Aligned.toTranscriptome.out.bam
rm ${OutputPath}/trim/${Sample}_RNA_P_1.fastq.gz
rm ${OutputPath}/trim/${Sample}_RNA_P_2.fastq.gz 

echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : RSEM Done \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log

wait


done < ${coldata}


wait
if [ ${dofeaturecount} = "y" ]
then

CA featurecount

wait 

echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : FeatureCounts \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log


    if [ ${dicision} = "UnStranded" ]
    then

        featureCounts -T 10 -s 0 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 2>> ${OutputPath}/log/${Sample}.log



    elif [ ${dicision} = "Stranded" ]
    then

        featureCounts -T 10 -s 1 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 2>> ${OutputPath}/log/${Sample}.log


    elif [ ${dicision} = "ReverselyStranded" ]
    then

        featureCounts -T 10 -s 2 -p \
        -a /users/data/reference/hg38/Homo_sapiens.GRCh38.104.gtf \
        -o ${OutputPath}/FeatureCounts/featurecounts.results.txt \
        ${OutputPath}/bam/*.STAR.Aligned.sortedByCoord.out.bam 2>> ${OutputPath}/log/${Sample}.log


    else
        echo "ERROR : undefined RSEM strandedness" |tee -a ${OutputPath}/log/${Sample}.log
        exit 1
    fi &&

    sed -i "s:${OutputPath}/bam/::" ${OutputPath}/FeatureCounts/featurecounts.results.txt
    sed -i "s:.STAR.Aligned.sortedByCoord.out.bam::" ${OutputPath}/featureCounts/featurecounts.results.txt
    cut -f1,7- ${OutputPath}/FeatureCounts/featurecounts.results.txt > ${OutputPath}/featureCounts/featurecounts.results.final.txt
elif [ ! ${dofeaturecount} = "n" ]
then
    echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
    date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
    echo -e "${Sample} : All Done w/o FeatureCounts \n\n" |tee -a ${OutputPath}/log/${Sample}.log
    echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
else
    echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
    date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
    echo -e "ERROR : n = do featurecounts, y = don't \n\n" |tee -a ${OutputPath}/log/${Sample}.log
    echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
    exit 1
fi &&

echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log
date +"%d-%m-%Y %T" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "${Sample} : All Done w/ FeatureCounts \n\n" |tee -a ${OutputPath}/log/${Sample}.log
echo -e "---------------------------------------------------- \n\n" |tee -a ${OutputPath}/log/${Sample}.log

#fin.