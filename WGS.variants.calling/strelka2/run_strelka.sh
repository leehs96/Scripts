#!/bin/bash
set -e

while IFS= read -r line || [ -n "$line" ] 
do
Line=$line

tumor=`echo ${Line} | awk '{print $1}'`
normal=`echo ${Line} | awk '{print $2}'`

COUNT_LOG_FILES=$(find . -maxdepth 1 -name "${tumor}}*" -type f | wc -l)

if [[ $COUNT_LOG_FILES -gt 0 ]]
then
rm /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}*
fi

if [ -d /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta ]
then
  rm -rf /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta
fi

mkdir /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta

echo -e "------------------------------ \n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
date +"%d-%m-%Y %T: ${tumor} manta START" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
echo -e "\n------------------------------ \n\n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log

configManta.py \
    --tumorBam /users/hslee/cancer/WGS_hg38/bam/${tumor}.s.md.br.bam \
    --normalBam /users/hslee/cancer/WGS_hg38/bam/${normal}.s.md.br.bam \
    --referenceFasta /users/hslee/ref/hg38/hg38_v0_Homo_sapiens_assembly38.fasta \
    --callRegions /users/hslee/ref/hg38/intervals/strelka2.hg38.callable.bed.gz \
    --runDir /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta &&

/users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta/runWorkflow.py \
    -m local \
    -j 20 &>> /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.manta.runWorkflow.out || { c=$?;echo "Error";exit $c; }

echo -e "------------------------------ \n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
date +"%d-%m-%Y %T: ${tumor} manta DONE" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
echo -e "\n------------------------------ \n\n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log

wait

if [ -d /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_strelka ]
then
  rm -rf /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_strelka
fi

mkdir /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_strelka

echo -e "------------------------------ \n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
date +"%d-%m-%Y %T: ${tumor} strelka START" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
echo -e "\n------------------------------ \n\n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log

configureStrelkaSomaticWorkflow.py \
    --reportEVSFeatures \
    --tumorBam /users/hslee/cancer/WGS_hg38/bam/${tumor}.s.md.br.bam \
    --normalBam /users/hslee/cancer/WGS_hg38/bam/${normal}.s.md.br.bam \
    --ref /users/hslee/ref/hg38/hg38_v0_Homo_sapiens_assembly38.fasta \
    --callRegions /users/hslee/ref/hg38/intervals/strelka2.hg38.callable.rm.cen.tel.gap.par.bed.gz \
    --indelCandidates /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_manta/results/variants/candidateSmallIndels.vcf.gz \
    --runDir /users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_strelka &&

/users/hslee/cancer/WGS_hg38/Illumina_pipe/${tumor}_strelka/runWorkflow.py \
    -m local \
    -j 20 &>> /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.strelka2.runWorkflow.out || { c=$?;echo "Error";exit $c; }

echo -e "------------------------------ \n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
date +"%d-%m-%Y %T: ${tumor} strelka done" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log
echo -e "\n------------------------------ \n\n" |tee -a /users/hslee/cancer/WGS_hg38/Illumina_pipe/log/${tumor}.log

done < /users/hslee/cancer/WGS_hg38/Illumina_pipe/name.txt