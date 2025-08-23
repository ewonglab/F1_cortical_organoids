### Generating mouse strain-specific genomes

Strain-specific variant call files were downloaded from the Mouse Genome Project repository - https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/ - including:
$$
\begin{array}{|c|c|}
\hline
\text{CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz} & \text{CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz} \\
\text{MOLF_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz} & \text{MOLF_EiJ.mgp.v5.snps.dbSNP142.vcf.gz} \\
\text{PWK_PhJ.mgp.v5.indels.dbSNP142.normed.vcf.gz} & \text{PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz} \\
\text{SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz} & \text{SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz} \\
\hline
\end{array}
$$
    
We used C57BL/6J mouse genome (GRCm38, release 68) as a reference genome. 
Downloaded from - https://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz 
Genome annotation - https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
