####################################################################################################
# Configuration

THREADS=16
PROJECT=/home/murilo/Documents/GIZ_Bolivia_2023/EjemploDENV
SAMPLES="SRR24053947"

####################################################################################################
# Activating environment and exporting configuration

conda deactivate
conda activate ngs
mamba env export > $PROJECT/env/ngs.yaml

####################################################################################################
# Pre-trimming quality

input=$PROJECT/raw/*.fastq.gz
output=$PROJECT/qc_before

fastqc -t $THREADS -o $output $input

####################################################################################################
# trimming

input=$PROJECT/raw
outputA=$PROJECT/trimmed
outputB=$PROJECT/qc_after

for sample in $SAMPLES
do
    R1=$input/${sample}_R1.fastq.gz
    R2=$input/${sample}_R2.fastq.gz

    trim_galore -j 8 -paired $R1 $R2 -o $outputA --fastqc_args "-t $THREADS --outdir $outputB/"

done

####################################################################################################
# Indexing

input=$PROJECT/index

bwa index $input/DENV_REF.fasta

####################################################################################################
# Mapping

index=$PROJECT/index/DENV_REF.fasta
input=$PROJECT/trimmed
output=$PROJECT/mapped

for sample in $SAMPLES
do
    R1=$input/${sample}_R1_val_1.fq.gz
    R2=$input/${sample}_R2_val_2.fq.gz
    BAM=$output/${sample}.bam

    bwa mem $index -t $THREADS $R1 $R2 | samtools view -bS -F 4 | samtools sort -o $BAM
    samtools index $BAM

done

####################################################################################################
# Control Coverage

input=$PROJECT/mapped
output=$PROJECT/coverage

for sample in $SAMPLES
do
    BAM=$input/${sample}.bam
    cov=$output/${sample}.samtools.coverage

    samtools coverage $BAM -m > $cov

done

####################################################################################################
# Consensus without primer removal

input=$PROJECT/mapped
output=$PROJECT/consensus

for sample in $SAMPLES
do
    BAM=$input/${sample}.bam
    con=$output/${sample}.notrim.fasta

    samtools consensus -a -m "bayesian" -o $con $BAM

done

####################################################################################################
# Consensus with iVar (primer clipping)

input=$PROJECT/mapped
primers=$PROJECT/primers/primers.bed
ref=$PROJECT/index/DENV_REF.fasta
output=$PROJECT/consensus

for sample in $SAMPLES
do
    BAM=$input/${sample}.bam
    tBAM=$input/${sample}.ptrim.bam

    con=$output/${sample}.trim.fasta

    ivar trim -e -i $BAM -b $primers | samtools sort -o $tBAM
    samtools mpileup -A -d 6000000 -B -Q 20 --reference $ref $tBAM | ivar consensus -p $con -m 10 -n N

done

####################################################################################################
