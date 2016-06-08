'''
TACOCAT: Trim And COnstruct Combined Assembly Transcriptome

This workflow takes raw Illumina paired-end fastq files and outputs a reference
transcriptome comprising 'primary' and 'alternate' sets of transcripts.

1. Quality trim reads with Trimmomatic.

2. Perform de novo assembly using Trinity.

3. Perform de novo assembly using Velvet and Oases under different k-mer sizes.

4. Merge all assemblies into a final reference assembly with EvidentialGene.
'''




__author__ = "Wade R. Roberts (http://wrroberts.github.io)"
__license__ = "GPL"




from os.path import join




FASTQ_DIR = './data/'




SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample, Samp[^/]+}_R1.fq'))




PATTERN_R1 = '{sample}_R1.fq'
PATTERN_R2 = '{sample}_R2.fq'




''' Quality trims Illumina paired-end reads taken from a fastq file '''
rule trimmomatic:
	input: 
		left=join(FASTQ_DIR, PATTERN_R1), 
		right=join(FASTQ_DIR, PATTERN_R2)
	output:
		paired1="trimmomatic/{sample}_R1_paired.fq",
		paired2="trimmomatic/{sample}_R2_paired.fq",
		unpaired1="trimmomatic/{sample}_R1_unpaired.fq",
		unpaired2="trimmomatic/{sample}_R2_unpaired.fq"
	benchmark:
		"benchmarks/{sample}.trimmomatic.benchmark.txt"
	threads: 8
	shell: 
		"java -jar /path/to/trimmomatic-0.36.jar PE -threads {threads} {input.left} {input.right} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} ILLUMINACLIP:/adapters/TruSeq2-PE.fa:2:30:10:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"



	
''' Combines all left reads from Trimmomatic output into a single fastq file '''
rule combine_trimmomatic_left:
	input:
		fq1P="trimmomatic/{sample}_R1_paired.fq",
		fq1U="trimmomatic/{sample}_R1_unpaired.fq"
	output:
		"trimmomatic/{sample}.reads.left.fq"
	shell:
		"cat {input.fq1P} {input.fq1U} > {output}"




''' Combines all right reads from Trimmomatic output into a single fastq file '''
rule combine_trimmomatic_right:
	input:
		fq2P="trimmomatic/{sample}_R2_paired.fq",
		fq2U="trimmomatic/{sample}_R2_unpaired.fq"
	output:
		"trimmomatic/{sample}.reads.right.fq"
	shell:
		"cat {input.fq2P} {input.fq2U} > {output}"
		
		
		

''' De novo Trinity assembly of trimmed paired-end reads using FR strand orientation and k-mer size 25 '''
rule trinity:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_trinity/trinity.fasta"
	benchmark:
		"benchmarks/{sample}.trinity.benchmark.txt"
	threads: 4
	shell:
		"perl /path/to/Trinity --seqType fq --max_memory 8G --left {input.left} --right {input.right} --SS_lib_type FR --CPU {threads} --full_cleanup --output {output}"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 25 '''
rule velvet_25:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/25"
	benchmark:
		"benchmarks/{sample}.velvet25.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 25 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 35 '''
rule velvet_35:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/35"
	benchmark:
		"benchmarks/{sample}.velvet35.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 35 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 45 '''
rule velvet_45:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/45"
	benchmark:
		"benchmarks/{sample}.velvet45.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 45 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 55 '''
rule velvet_55:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/55"
	benchmark:
		"benchmarks/{sample}.velvet55.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 55 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 65 '''
rule velvet_65:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/65"
	benchmark:
		"benchmarks/{sample}.velvet65.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 65 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo Velvet assembly of trimmed paired-end reads with a k-mer size 75 '''
rule velvet_75:
	input:
		left="trimmomatic/{sample}.reads.left.fq",
		right="trimmomatic/{sample}.reads.right.fq"
	output:
		"{sample}_oases_out/75"
	benchmark:
		"benchmarks/{sample}.velvet75.benchmark.txt"
	shell:
		"/path/to/velveth {output}/ 75 -fastq -separate -shortPaired {input.left} {input.right} && /path/to/velvetg {output}/ -read_trkg yes -ins_length 300 -min_contig_lgth 200 -unused_reads yes -clean yes"




''' De novo assembly with Oases using the Velvet output at k-mer size 25 '''
rule oases_25:
	input:
		"{sample}_oases_out/25"
	output:
		"{sample}_oases_out/25/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases25.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' De novo assembly with Oases using the Velvet output at k-mer size 35 '''
rule oases_35:
	input:
		"{sample}_oases_out/35"
	output:
		"{sample}_oases_out/35/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases35.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' De novo assembly with Oases using the Velvet output at k-mer size 45 '''
rule oases_45:
	input:
		"{sample}_oases_out/45"
	output:
		"{sample}_oases_out/45/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases45.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' De novo assembly with Oases using the Velvet output at k-mer size 55 '''
rule oases_55:
	input:
		"{sample}_oases_out/55"
	output:
		"{sample}_oases_out/55/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases55.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' De novo assembly with Oases using the Velvet output at k-mer size 65 '''
rule oases_65:
	input:
		"{sample}_oases_out/65"
	output:
		"{sample}_oases_out/65/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases65.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' De novo assembly with Oases using the Velvet output at k-mer size 75 '''
rule oases_75:
	input:
		"{sample}_oases_out/75"
	output:
		"{sample}_oases_out/75/transcripts.fa"
	benchmark:
		"benchmarks/{sample}.oases75.benchmark.txt"
	shell:
		"/path/to/oases {input}/ -ins_length2 300 -min_trans_lgth 200"




''' Concatenate all Oases assemblies into a single fasta file '''
rule combine_oases_assemblies:
	input:
		f25="{sample}_oases_out/25/transcripts.fa",
		f35="{sample}_oases_out/35/transcripts.fa",
		f45="{sample}_oases_out/45/transcripts.fa",
		f55="{sample}_oases_out/55/transcripts.fa",
		f65="{sample}_oases_out/65/transcripts.fa",
		f75="{sample}_oases_out/75/transcripts.fa"
	output:
		"{sample}_oases_out/all.transcripts.oases.fasta"
	shell:
		"cat {input.f25} {input.f35} {input.f45} {input.f55} {input.f65} {input.f75} > {output}"




''' Reformat Oases headers to a standard format needed by EviGene '''
rule change_oases_format:
	input:
		"{sample}_oases_out/all.transcripts.oases.fasta"
	output:
		"{sample}_oases_out/all.transcripts.oases.trformat.fasta"
	shell:
		"perl /path/to/trformat.pl -input {input} -output {output} -format=Oases -MINTR=200 -prefix Oases"




''' Reformat Trinity headers to a standard format needed by EviGene '''
rule change_trinity_format:
	input:
		"{sample}_trinity/trinity.fasta"
	output:
		"{sample}_trinity/trinity.trformat.fasta"
	shell:
		"perl /path/to/trformat.pl -input {input}.Trinity.fasta -output {output} -format=Trinity -MINTR=200 -prefix Trinity"




''' Merge Trinity and Oases assemblies into primary and alternate transcript sets using EviGene '''
rule evigene:
	input:
		oases="{sample}_oases_out/all.transcripts.oases.trformat.fasta",
		trinity="{sample}_trinity/trinity.trformat.fasta"
	output:
		"EviGene/{sample}.evigene.fasta"
	benchmark:
		"benchmarks/{sample}.evigene.benchmark.txt"
	threads: 8
	shell:
		"cat {input.oases} {input.trinity} > {output} && mkdir EviGene/results && cd EviGene/results && perl /path/to/tr2aacds.pl -mrnaseq {output} -MINCDS=75 -NCPU {threads} -MAXMEM=60000 -logfile -tidyup"



