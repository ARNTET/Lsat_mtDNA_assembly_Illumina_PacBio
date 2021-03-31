#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000

# Modules needed ###############################################################
module load cutadapt
module load NOVOPlasty
module load bwa
module load samtools
module load canu

# /!\ To fill /!\ ##############################################################
samples=()
lengthS=${#samples[@]}
adapter_FWD=
adapter_REV=
config_files=()

# Step one #######################################################################
a=0
while (($a<$lengthS)); do
  #A# Trimming
  cutadapt \
  -a $adapter_FWD \
  -A $adapter_REV \
  -o ${samples[a]}".R1.fq.gz" \
  -p ${samples[a]}".R2.fq.gz" \
  ${samples[a]}"_R1.fastq.gz" \
  ${samples[a]}"_R2.fastq.gz"
  perl NOVOPlasty3.2.pl \
  -c ${config_files[a]}
  samtools faidx b_assemblyN/${samples[a]}".novoplasty"/contigs.fasta
  wc -l b_assemblyN/${samples[a]}".novoplasty"/contigs.fasta.fai
  let "a=a+1"
done
echo Now contigs need to be blasted to identify mitochondrial coverage

# /!\ To fill /!\ ##############################################################
a=0
while (($a<$lengthS)); do
echo What is the upper limit to mitochondrial coverage in ${samples[a]}
read ${samples[a]}"_Upper"
echo What is the lower limit to mitochondrial coverage in ${samples[a]}
read ${samples[a]}"_Lower"
let "a=a+1"
done
echo Thank you for these informations

# Step two #######################################################################
a=0
while (($a<$lengthS)); do
  #C# Filtering
  mkdir c_filtering
  cat b_assemblyN/${samples[a]}".novoplasty"/contigs.fasta.fai \
  | sed -e $'s/_/\t/g' \
  | sed 's/\./,/g' \
  | awk '$6 < '${samples[a]}"_Upper"' && $6 > '${samples[a]}"_Lower"'' \
  | cut -f 1-6 \
  | sed -e $'s/\t/_/g' \
  | sed 's/,/\./g' > c_filtering/${samples[a]}"list_contigs.txt"
  xargs samtools faidx \
  b_assemblyN/${samples[a]}".novoplasty"/contigs.fast \
  < c_filtering/${samples[a]}"list_contigs.txt" \
  > c_filtering/${samples[a]}"contigs_flt.fasta"
  #D# Mapping of Pacbio reads on Illumina contigs
  mkdir d_Pacbio_mapping
  bwa index c_filtering/${samples[a]}"contigs_flt.fasta"
  bwa mem -t 16 \
  c_filtering/${samples[a]}"contigs_flt.fasta" \
  ${#samples[@]}".pacbio.fastq" \
  > d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.sam"
  samtools view -b \
  d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.sam" \
  > d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.bam"
  samtools fastq -F4 \
  d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.bam" \
  > d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.fastq"
  #D# Assembly of contigs from Pacbio mapped reads
  mkdir e_assemblyC
  canu \
  -p lsat \
  -d e_assemblyC \
  genomeSize=0.3m \
  -pacbio-raw d_Pacbio_mapping/${#samples[@]}".mapped.pacbio.fastq"
  let "a=a+1"
done
