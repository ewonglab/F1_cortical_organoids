workdir: "/projects/F1_organoids"
configfile: "github/cluster.yml"

strain_names = ["CAST_EiJ", "MOLF_EiJ", "PWK_PhJ", "SPRET_EiJ"]
strain_list = ["CAST", "MOLF", "PWK", "SPRET"]

from datetime import date
today = date.today()
dt = today.strftime("%Y%m%d")


samples_table = pd.read_csv("github/rna_samples.csv").set_index("sample", drop=False)

SAMPLES = list(samples_table['sample'])

ASSAYS = list(samples_table['assay'])

STRAINS = list(samples_table['strain'])

rule all:
	input:
		expand(["counts/GRCm38/{assay}/{strain}/{sample}/{assay}_{sample}_gene_counts.txt"], zip, assay = ASSAYS, sample = SAMPLES, strain = STRAINS)
		
		
####################
# FASTQC RNA files #
####################


def get_fq1(wildcards):
	return glob.glob('fastq/'+ wildcards.assay + "/" + wildcards.sample + "/" + '*_' + wildcards.run + '*R1_001.fastq.gz') 

def get_fq2(wildcards):
	return glob.glob('fastq/'+ wildcards.assay + "/" + wildcards.sample + "/" + '*_' + wildcards.run + '*R2_001.fastq.gz') 

rule fastqc_raw_rna:
	input:
		get_fq1,
		get_fq2
	output:
		zip  = "fastq/qc/{assay}/{sample}/{sample}_fastqc.zip",
		html = "fastq/qc/{assay}/{sample}/{sample}_fastqc.html",
	params:
		dir = "fastq/qc/{assay}/{sample}/}",
		tools = "/tools/FastQC",
		tmp = "tmp",
		log = "fastq/qc/{assay}/{sample}/fastqc_{sample}.log",
	threads: 4
	shell:
		"""
		{params.tools}/fastqc -o {params.dir} -d {params.tmp} {input.fq1} {input.fq2} 2> {params.log}
		"""

###########################
# Trimming bulk RNA fastq #
###########################

rule trim_fastq:
	input:
		fq1= "fastq/{assay}/{sample}/R1.fastq.gz",
		fq2= "fastq/{assay}/{sample}/R2.fastq.gz",
	output:
		fq1_trim = "fastq/{assay}/{sample}/R1.trimmom.fastq.gz",
		fq1_trim_unp = "fastq/{assay}/{sample}/R1.trimmom.unp.fastq.gz",
		fq2_trim = "fastq/{assay}/{sample}/R2.trimmom.fastq.gz",
		fq2_trim_unp = "fastq/{assay}/{sample}/R2.trimmom.unp.fastq.gz"
	params:
		bin = "/tools/Trimmomatic-0.39",
		log = "logs/trimmomatic_{assay}_{sample}.log"
	shell:
		"""
		java -jar {params.bin}/trimmomatic-0.39.jar PE -threads 4 -trimlog {params.log} \
		{input} \
		{output} \
		ILLUMINACLIP:{params.bin}/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:30
		"""


####################
# Alignng bulk RNA #
####################

rule star_align_rna:
	input:
		fq1_trim = "fastq/{assay}/{sample}/R1.trimmom.fastq.gz",
		fq2_trim = "fastq/{assay}/{sample}/R2.trimmom.fastq.gz",
		genome_dir = "genomes/combined/{strain}/star_index_GRCm38"
	output:
		align_bam = directory("align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_")
	params:
		software = "/g/data/zk16/software/STAR-2.7.6a/bin/Linux_x86_64",
		logs = "logs/{assay}_{sample}_combined_{strain}_star.log"
	threads: 12
	shell:
		"""
		{params.software}/STAR --runThreadN {threads} --runMode alignReads --genomeDir {input.genome_dir} --readFilesIn {input.fq1_trim} {input.fq2_trim} --readFilesCommand zcat \
		--outFileNamePrefix {output} --outSAMtype BAM SortedByCoordinate --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMatchNmin 30 \
		--alignIntronMin 20 --alignIntronMax 20000 --alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM --outFilterMismatchNoverReadLmax 0.001 --scoreDelOpen -1000 --scoreInsOpen -1000 > {params.logs} 2>&1
		"""
		

################################
# Subset uniquely mapped reads #
################################

rule bam_index:
	input:
		align_bam = directory("align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_"),
	output:
		bam_nomismatch = "align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_nomismatch.sortedByCoord.bam"
	shell:
		"""
		samtools view -h {input} | grep -P '^@|NM:i:0\b' | samtools view -Sb - > {output}
		samtools index {output}
		"""
		
		
########################################################
# Remove reads mapping to black-listed variant regions #
########################################################		

rule bam_blregion:
	input:
		bam_nomismatch = "align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_nomismatch.sortedByCoord.bam",
		bl_region = "/variants/black_regions_mm10/{strain}_bl_region_forolap.bed"
	output:
		bam_bl_list_rm = "align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_noNM_bllistrm.bam",
		bam_bl_list = "align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_noNM_bllist.bam"
	threads: 12
	params:
		output_dir = "align/GRCm38/{assay}/{strain}/{sample}/"
	shell:
		"""
		mkdir -p {params.output_dir}
		samtools view {input.bam_nomismatch} -@ {threads} -h -b -o {output.bam_bl_list} -L {input.bl_region} -U {output.bam_bl_list_rm}
		samtools index {output.bam_bl_list_rm}
		"""


##########################
# Allele-specific counts #
##########################		


rule allelic_counts:
	input:
		bam_bl_list_rm = "align/GRCm38/{assay}/{strain}/{sample}/combined_{assay}_{sample}_noNM_bllistrm.bam",
	output:
		counts = "counts/GRCm38/{assay}/{strain}/{sample}/{assay}_{sample}_gene_counts.txt",
	params:
		strain_id = "{strain}",
		ref_gtf = "genomes/reference/Mus_musculus.GRCm38.102.subset.gtf",
		pseudo_gtf = "genomes/SPRET/SPRET_pseudo_v102_v5_nomask.gtf",
		feature = "exon", #select from either gene or exon 
		Rlog = "logs/bulkRNA_day12_organoids_{sample}_allelic_counts.log",			
	shell:
		"""
		mkdir -p {params.dir}
		Rscript R/gene_counts_mm10.R {params.ref_gtf} {params.pseudo_gtf} {params.strain_id} {params.feature} {input.bam_bl_list_rm} {output.counts} >& {params.Rlog}
		"""


		