#!/bin/bash
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name="Illumina_Pa_Assembly"
#SBATCH --output=JGAF-Illumina_Assembly%j.out

# Modules needed ###############################################################
module load NOVOPlasty
module load fastq-join
module load samtools
module load canu

# Sample #######################################################################

#A# ASSEMBLY of contigs from Illumina reads
mkdir a_Illumina_assembly
perl NOVOPlasty3.2.pl -c config.txt

#B# Manual identification of mitochondrial contigs
mkdir b_Illumina_mt
bwa index mt_Illumina_contigs.fasta

#C# Mapping of Pacbio reads on Illumina contigs
mkdir c_Pacbio_mapping
bwa mem -t 16 \
b_Illumina_mt/t_Illumina_contigs.fasta \
pacbio.fastq \
> c_Pacbio_mapping/mapped.pacbio.sam
samtools view -b \
c_Pacbio_mapping/mapped.pacbio.sam \
> c_Pacbio_mapping/mapped.pacbio.bam
samtools fastq -F4 \
c_Pacbio_mapping/mapped.pacbio.bam \
> c_Pacbio_mapping/mapped.pacbio.fastq

#D# Assembly of contigs from Pacbio mapped reads
mkdir d_assemblyC
canu \
 -p lsat \
 -d d_assemblyC \
 genomeSize=0.3m \
 -pacbio-raw mapped.pacbio.fastq
