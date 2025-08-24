workdir: "/projects/F1_organoids"
configfile: "github/cluster.yml"

strain_names = ["CAST_EiJ", "MOLF_EiJ", "PWK_PhJ", "SPRET_EiJ"]
strain_list = ["CAST", "MOLF", "PWK", "SPRET"]

from datetime import date
today = date.today()
dt = today.strftime("%Y%m%d")


rule all:
	input:
		filter_vars = expand(["variants/{strain_id}.mgp.v5.snp.indels.sorted.PASS.vcf.gz"], strain_id = strain_names)

		
###########################
# preprocessing vcf files #
###########################

rule merge_vcf:
	input:
		all_snps = "variants/{strain_id}.mgp.v5.snps.dbSNP142.vcf.gz",
		all_indels = "variants/{strain_id}.mgp.v5.indels.dbSNP142.normed.vcf.gz",
	output:
		merged_vars = "variants/{strain_id}.mgp.v5.snp.indels.dbSNP142.vcf.gz", 
	threads: 12
	shell:
		"""
		bcftools concat -n --threads {threads} -o {output.merged_vars} {input.all_snps} {input.all_indels} 
		"""

rule filter_vcf:
	input:
		merged_vars = "variants/{strain_id}.mgp.v5.snp.indels.dbSNP142.vcf.gz", 
	output:
		filter_vars = "variants/{strain_id}.mgp.v5.snp.indels.sorted.PASS.vcf.gz", 
	shell:
		"""
		bcftools view -f PASS -s {strain_id} {input.merged_vars} | \
		bcftools sort -Oz - > {output.filter_vars}
		"""
