# TACOCAT version 0.0.1 - A pipeline for de novo transcriptome assembly.

Available from:

http://github.com/wrroberts/TACOCAT


**TACOCAT** - **T**rim **A**nd **CO**nstruct **C**ombined **A**ssembly **T**ranscriptome

TACOCAT was built to streamline the de novo transcriptome assembly process and create a robust reference for downstream analysis. This pipeline starts with raw paired-end Illumina reads in fastq format, quality trims the reads, assembles the reads using multiple assemblers and different k-mer sizes, and outputs a reference classified into primary and alternate sets of transcripts.


Usage
=====
TACOCAT runs as a single command that takes as input paired-end Illumina fastq files, two per species, outputs a number of directories containing intermediate files from each step of the pipeline, and provides the final reference transcriptome classified into primary and alternate transcripts. Classification is based on CDS-dna local alignment identity. Perfect fragment CDS are dropped, those with some CDS base differences are kept, with the longest CDS as primary transcript. Alternates are any alternately spliced transcripts.

INPUT - two fastq files. One file for forward reads and one file for reverse reads. If you have multiple forward and reverse read files from  replicate conditions or sequencing runs, combine these together before using the pipeline.

If you want to run the pipeline using multiple processors in parallel, modify `config.yaml` to specify the number of processors you wish to use.

Installing Programs and Dependencies
====================================
TACOCAT pipeline is written to run on linux using the Snakemake workflow management system and requires the following programs and dependencies to be installed and in the system path:

1. Python 3.x 

2. Snakemake

3. Trimmomatic v.0.36

4. Trinity

5. Velvet v.1.2.10

6. Oases v.0.2.08

7. EvidentialGene

8. bowtie

9. fastanrdb

10. cd-hit, cd-hit-est

11. BLAST+

Brief instructions are given below although users may wish to refer to the installation notes provided with these packages for more detailed instructions. 

**python**

https://www.python.org.

**Snakemake**

Sourcecode and documentation can be found at https://bitbucket.org/snakemake/snakemake/src.

**Trimmomatic**

Binary sourcecode and documentation can be found at https://www.usadellalab.org/cms/?page=trimmomatic. Trimmomatic requires Java to access and run jar files. 

**Trinity**

Sourcecode and documentation can be found at https://github.com/trinityrnaseq/trinityrnaseq.

**Velvet**

Sourcecode and manual can be found at https://www.ebi.ac.uk/~zerbino/velvet/.

**Oases**

Sourcecode and manual can be found at https://www.ebi.ac.uk/~zerbino/oases/.

**EvidentialGene**

Sourcecode can be found at https://sourceforge.net/projects/evidentialgene/. Additional documentation about the program can be found at https://arthropods.eugenes.org/EvidentialGene/. TACOCAT uses `trformat.pl` (/scripts/rnaseq/) and `tr2aacds.pl` (/scripts/prot/) scripts.

**bowtie**

Sourcecode and documentation can be found at https://bowtie-bio.sourceforge.net/index.shtml.

**fastanrdb**

Part of the Exonerate package. Sourcecode and documentation can be found at https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate.

**cd-hit, cd-hit-est**

Sourcecode and documentation can be found at https://weizhongli-lab.org/cd-hit/.

**BLAST+**

Executables are found here ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST (instructions are currently at http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/ and in more detail in the BLAST+ user manual). As websites can change, an alternative is to search online for "install BLAST+".

1. Instructions are provided for installing BLAST+ on various flavors of linux on the 'Standalone BLAST setup for Unix' page of the BLAST+ Help manual currently at http://www.ncbi.nlm.nih.gov/books/NBK1762/.
2. Follow the instructions under "Configuration" in the BLAST+ help manual to add BLAST+ to the PATH environment variable.

Setting up and running TACOCAT
===============================
Once the required programs and dependencies have been installed, TACOCAT can be setup and run on the small example data included in the package as follows:

1. Save TACOCAT_v0.0.1.tar.gz
2. Open a terminal and cd to the directory where TACOCAT_v.0.0.1.tar.gz was saved
3. `tar xzf TACOCAT_v.0.0.1`
4. `cd TACOCAT_v.0.0.1`
5. `snakemake EviGene/{sample}.evigene.fasta`

The command for running TACOCAT on the example dataset is:

`snakemake EviGene/{A,B}.evigene.fasta `

