workdir: "/projects/F1_organoids"
configfile: "github/cluster.yml"

strain_names = ["CAST_EiJ", "MOLF_EiJ", "PWK_PhJ", "SPRET_EiJ"]
strain_list = ["CAST", "MOLF", "PWK", "SPRET"]

from datetime import date
today = date.today()
dt = today.strftime("%Y%m%d")


rule all:
    input:
        genome_dir = expand("genomes/combined/{strain}/star_index_GRCm39_PASS_nomask", strain = strain_list)

		
#########################
# make pseudo-reference #
#########################

rule pseudo_reference:
	input:
		filter_vars = "variants/{strain_id}.mgp.v5.snp.indels.sorted.PASS.vcf.gz", 
		fasta_ref = "genomes/reference/GRCm38_v68.genome.sorted.fa",	
	output:
		pseudo_ref = "genomes/{strain}/{strain}_GRCm38_v68_pseudo_PASS.fa",
		chain = "genomes/{strain}/GRCm38_v68_to_{strain}_v5_PASS.chain",
	shell:
		"""
		bcftools consensus {input.filter_vars} --sample {strain_id} --fasta-ref {input.fasta_ref} -c {output.chain} -o {output.pseudo_ref}
		"""
		
rule pseudo_gtf:
	input:
		ref_gtf = "genomes/reference/Mus_musculus.GRCm38.102.subset.gtf",
		chain = "genomes/{strain}/GRCm38_v68_to_{strain}_v5_PASS.chain",
	output:
		pseudo_gtf = "genomes/{strain}/{strain}_pseudo_GRCm38_v102_PASS.gtf",
		pseudo_gtf_unmp = "genomes/{strain}/{strain}_pseudo_GRCm38_v102_PASS.gtf.unmapped",
	params:
		tools = "tools/utils",
	shell:
		"""
		{params.tools}/liftOver -gff {input.ref_gtf} {input.chain} {output.pseudo_gtf} {output.pseudo_gtf_unmp}
		"""

rule rename_gtf:
	input:
		pseudo_gtf = "genomes/{strain}/{strain}_pseudo_GRCm38_v102_PASS.gtf",
		pseudo_ref = "genomes/{strain}/{strain}_GRCm38_v68_pseudo_PASS.fa",
	output:
		renamed_pseudo_ref = "genomes/{strain}/{strain}_GRCm38_v68_pseudo_PASS_renamed.fa",
		renamed_pseudo_gtf = "genomes/{strain}/{strain}_pseudo_GRCm38_v102_PASS_renamed.gtf",
	shell:
		"""
		sed 's/^>/>pseudo/' {input.pseudo_ref} > {output.renamed_pseudo_ref}
		awk '{ $1="pseudo" $1; } 1' {input.pseudo_gtf} > {output.renamed_pseudo_gtf}
		"""
		
rule comb_reference:
	input:
		ref_gtf = "genomes/reference/Mus_musculus.GRCm38.102.subset.gtf",
		fasta_ref = "genomes/reference/GRCm38_v68.genome.sorted.fa",	
		renamed_pseudo_ref = "genomes/{strain}/{strain}_GRCm38_v68_pseudo_PASS_renamed.fa",
		renamed_pseudo_gtf = "genomes/{strain}/{strain}_pseudo_GRCm38_v102_PASS_renamed.gtf",
	output:
		final_gtf = "genomes/combined/{strain}/GRCm38_v102_{strain}_PASS.gtf",
		final_ref = "genomes/combined/{strain}/GRCm38_v168_{strain}_PASS.fa"
	shell:
		"""
		cat {input.ref_gtf} {input.renamed_pseudo_gtf} > {output.final_gtf}
		cat {input.fasta_ref} {input.renamed_pseudo_ref} > {output.final_ref}
		"""
		
#####################
# create STAR index #
#####################
		
rule star_index:
	input:
		final_gtf = "genomes/combined/{strain}/GRCm38_v102_{strain}_PASS.gtf",
		final_ref = "genomes/combined/{strain}/GRCm38_v168_{strain}_PASS.fa"
	output:
		genome_dir = "genomes/combined/{strain}/star_index_GRCm38"
	params:
		software = "software/STAR-2.7.6a/bin/Linux_x86_64",
		logs = "logs/genomeindex_GRCm38_combined_{strain}.log"
	threads: 16
	shell:
		"""
		{params.software}/STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.genome_dir} --genomeFastaFiles {input.final_ref} \
		--sjdbGTFfile {input.final_gtf} --sjdbOverhang 100 --genomeSuffixLengthMax 300 --limitGenomeGenerateRAM 64000000000 > {params.logs} 2>&1
		"""
		
		