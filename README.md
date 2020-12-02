Mitochondrial genome assembly combining Pacbio SMRT and Illumina MiSeq technologies

# De novo assembly of plant mitochondrial DNA using Pacbio SMRT and Illumina MiSeq technologies

## Materials

### Raw reads fastq files are available on ebi.ac.uk/ena under the accession numbers:

- Lactuca sativa:
  - Pacbio: ERS5266954
  - Illumina Miseq (2 x 151 bp): ERS5266978

## Methods

### Illumina sequencing
Library preparation and sequencing reactions were performed by Mr. A. Alioua and Ms. S. Koechler from the IBMP “gene expression analysis” platform. Libraries were prepared with the Nextera XT DNA Library Prep kit (Illumina, San Diego, California, USA) and their size checked on the 2100 Bioanalyzer (Agilent, Santa Clara, California, USA). The libraries were then sequenced with paired-end 2 x 151 reads on a MiSeq System (Illumina). The raw Illumina reads were deposited in GenBank under accession numbers ERS5266978 and ERS5267070 to ERS5267078 (www.ebi.ac.uk/ena).

### Pacbio sequencing
L. sativa var. capitata L. nidus tenerrima Pacbio sequencing was performed by Keygene (Wageningen, The Netherlands) on a Pacbio Sequel II instrument (Pacific Biosciences, Menlo Park, California, USA). The sequencing generated 19.9 Gb from 3.2 million reads (average read size 3051 pb).

### Mitochondrial assembly

Contigs from Illumina MiSeq paired-end (2 x 151 bases) were built with the organellar de novo assembler Novoplasty (k-mer 22). This led to 48 contigs, among them 22 were confirmed as mtDNA contigs by BLAST. Repeated regions were identified as the overlap between 4 different contigs and particularly three large ones: : R01 (34 696 bp), R02 (10 430 bp) and R03 (3 552 bp).

Raw Pacbio data delivered by Keygene were mapped on the 21 mtDNA contigs from Illumina with bwa-mem (see JGAF-Pacbio_Illumina_mtDNA_assembly.sh). Mapped Pacbio reads were then assembled using Canu to obtain 4 contigs. Among them 3 large ones were identified as mitochondrial by BLAST: 177 077 bp (PB01), 90 502 bp (PB02), and 73 475 bp (PB03).

The start of the PB01 sequence overlaps its end for 5 225 bp. Consequently, we circularized it into a circular chromosome (PBc01) of 171 852 bp. Furthermore, we found that the extremities of contigs PB02 and PB03 mapped to sequences internal to PBc01, corresponding to repeated regions identified among Novoplasty contigs. So we circularized PB02 and PB03 using repeats R01 and R02, respectively, into circular chromosomes PBc02 (112 968 bp) and PBc03 (78 504 bp). The circular chromosomes were then assembled through recombination-like events involving repeated sequences R01 and R02 into a single circular molecule of 363 324 bp.
