SAMPLES, = glob_wildcards('{sample}_R1_001.fastq.gz')

rule all:
    input:
        expand("filtered_reads/{sample}_nohost_R1.fastq", sample=SAMPLES),
        expand("filtered_reads/{sample}_nohost_R2.fastq", sample=SAMPLES),
        expand("kaiju_outputs/{sample}_kaiju.output", sample=SAMPLES),
        expand("kaiju_summaries/{sample}_table.tsv", sample=SAMPLES),
        expand("metaspades_output/{sample}/contigs.fasta", sample=SAMPLES),
       # expand("metaquast_output/{sample}/report.html", sample=SAMPLES),
        expand("assembly_database/{sample}.1.bt2", sample=SAMPLES),
        expand("mapped_to_assembly/{sample}_sorted.bam", sample=SAMPLES),
        expand("mapped_to_assembly/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("assembly_stats/{sample}_idxstats.txt", sample=SAMPLES),
       # expand("assembly_stats/{sample}_counts.txt", sample=SAMPLES),
        expand("centrifuge_reports/{sample}_centrifuge_report.tsv", sample=SAMPLES),
        expand("centrifuge_classifications/{sample}_classification.txt", sample=SAMPLES),
        expand("filtered_ids/{sample}_filtered_read_ids.txt", sample=SAMPLES),
	expand("filtered_reads_cyb/{sample}_host_cyb_R1.fastq", sample=SAMPLES),
	expand("filtered_reads_cyb/{sample}_host_cyb_R2.fastq", sample=SAMPLES),
	expand("metaspades_output_cyb/{sample}/contigs.fasta", sample=SAMPLES),
	expand("filtered_reads_bam/{sample}.bam.bai", sample=SAMPLES),
	expand("krona_charts/{sample}_krona.html", sample=SAMPLES),
	expand("extracted_reads/{sample}_R1_extracted_reads.fastq", sample=SAMPLES),
        expand("extracted_reads/{sample}_R2_extracted_reads.fastq", sample=SAMPLES)
# Since the genome is manually downloaded, the download rule is removed.

# Decompression step is removed if your genome is already in .fna (uncompressed) format.

# Build the Bowtie2 database
rule build_bowtie2_db:
    input:
    	"genome/host.fna"
    output:
    	"genome/host.1.bt2"
    shell:
    	"bowtie2-build {input} genome/host"

# Map reads against the host database and filter host reads
rule map_reads:
    input:
    	db="genome/host.1.bt2",
        R1="{sample}_R1_001.fastq.gz",
        R2="{sample}_R2_001.fastq.gz"
    output:
    	sam="mapped_reads/{sample}_host_mapped.sam"
    params:
        db_prefix="genome/host"
    shell:
    	"bowtie2 --very-sensitive-local -p 8 -x {params.db_prefix} -1 {input.R1} -2 {input.R2} -S {output}"

# Convert SAM to BAM
rule sam_to_bam:
    input:
    	"mapped_reads/{sample}_host_mapped.sam"
    output:
    	"mapped_reads/{sample}_host_mapped.bam"
    shell:
    	"samtools view -bS {input} > {output}"

# Extract unmapped reads
rule extract_unmapped_reads:
    input:
    	"mapped_reads/{sample}_host_mapped.bam"
    output:
    	temp("filtered_reads/{sample}_unmapped.bam")
    shell:
    	"samtools view -b -f 12 -F 256 {input} > {output}"

# Sort unmapped reads by name
rule sort_reads:
    input:
    	"filtered_reads/{sample}_unmapped.bam"
    output:
    	"filtered_reads/{sample}_unmapped_sorted.bam"
    shell:
    	"samtools sort -n {input} -o {output}"

# Convert sorted BAM to FastQ
rule bam_to_fastq:
    input:
    	"filtered_reads/{sample}_unmapped_sorted.bam"
    output:
    	R1="filtered_reads/{sample}_nohost_R1.fastq",
        R2="filtered_reads/{sample}_nohost_R2.fastq"
    shell:
    	"samtools fastq -1 {output.R1} -2 {output.R2} {input}"

#Kaiju
rule run_kaiju:
    input:
    	R1="filtered_reads/{sample}_nohost_R1.fastq",
        R2="filtered_reads/{sample}_nohost_R2.fastq"
    output:
    	"kaiju_outputs/{sample}_kaiju.output"
    params:
    	kaiju_db_nodes = "kaijudb/nodes.dmp",
        kaiju_db_fmi = "kaijudb/viruses/kaiju_db_viruses.fmi",
        threads =24,  # Adjust as needed
        e_value = "1e-05"
    shell:
    	"kaiju -t {params.kaiju_db_nodes} -f {params.kaiju_db_fmi} -i {input.R1} -j {input.R2} -o {output} -z {params.threads} -E {params.e_value} -v"

#Kaiju summary table
rule kaiju2table_summary:
    input:
    	"kaiju_outputs/{sample}_kaiju.output"
    output:
    	"kaiju_summaries/{sample}_table.tsv"
    params:
    	kaiju_db_nodes = "kaijudb/nodes.dmp",
        kaiju_db_names = "kaijudb/names.dmp",
        taxonomic_level = "species",  # Can be adjusted as needed
        taxonomic_levels = "superkingdom,phylum,class,order,family,genus,species"  # Adjust as needed
    shell:
    	"kaiju2table -t {params.kaiju_db_nodes} -n {params.kaiju_db_names} -r {params.taxonomic_level} -o {output} {input} -l {params.taxonomic_levels}"

#metaspades
rule run_metaspades:
    input:
        R1="filtered_reads/{sample}_nohost_R1.fastq",
        R2="filtered_reads/{sample}_nohost_R2.fastq"
    output:
        contigs="metaspades_output/{sample}/contigs.fasta",
        scaffolds="metaspades_output/{sample}/scaffolds.fasta"
    params:
        spades_path="spades.py",  # Update with your path to SPAdes
        threads=48,  # Adjust as needed
        memory=132,  # Adjust as needed, in GB
        outdir=lambda wc: "metaspades_output/" + wc.sample  # Dynamically generate output directory path
    shell:
        """
        {params.spades_path} --meta -1 {input.R1} -2 {input.R2} -k 55,75,95 \
        -t {params.threads} -m {params.memory} -o {params.outdir}
        """

#metaquast
rule evaluate_assembly:
    input:
        contigs="metaspades_output/{sample}/contigs.fasta"
    output:
        report="metaquast_output/{sample}/report.html"
    params:
        metaquast_path="metaquast.py",  # Update with your path to metaQUAST
        min_contig_length=10,  # Update as needed
        threads=24,  # Adjust as needed
        outdir=lambda wc: "metaquast_output/" + wc.sample  # Dynamically generate output directory path
    shell:
        """
        {params.metaquast_path} -o {params.outdir} \
        -m {params.min_contig_length} -t {params.threads} {input.contigs}
        """
rule run_centrifuge:
    input:
        R1="filtered_reads/{sample}_nohost_R1.fastq",
        R2="filtered_reads/{sample}_nohost_R2.fastq"
    output:
        report="centrifuge_reports/{sample}_centrifuge_report.tsv",
        classification="centrifuge_classifications/{sample}_classification.txt"
    params:
        centrifuge_index="centrifuge/index/p_compressed+h+v",  # Update with your path to Centrifuge index
        centrifuge_path="centrifuge",  # Update with your path to Centrifuge executable
        threads=34  # Adjust as needed
    shell:
        """
        {params.centrifuge_path} -x {params.centrifuge_index} \
        -1 {input.R1} -2 {input.R2} \
        -S {output.classification} \
        --report-file {output.report} \
        -p {params.threads} --min-hitlen 15
        """

rule build_assembly_database:
    input:
        "metaspades_output/{sample}/contigs.fasta"
    output:
        "assembly_database/{sample}.1.bt2"
    params:
        seed=4
    shell:
        "bowtie2-build --seed {params.seed} {input} assembly_database/{wildcards.sample}"
rule map_reads_to_assembly:
    input:
        db="assembly_database/{sample}.1.bt2",
        R1="filtered_reads/{sample}_nohost_R1.fastq",
        R2="filtered_reads/{sample}_nohost_R2.fastq"
    output:
        "mapped_to_assembly/{sample}.sam"
    params:
        threads=24,
        seed=4,
        db_prefix="assembly_database/{sample}"
    shell:
        "bowtie2 --sensitive-local -p {params.threads} --seed {params.seed} -x {params.db_prefix} -1 {input.R1} -2 {input.R2} -S {output}"
rule index_assembly:
    input:
        "metaspades_output/{sample}/contigs.fasta"
    output:
        "metaspades_output/{sample}/contigs.fasta.fai"
    shell:
        "samtools faidx {input}"
rule sam2_to_bam2:
    input:
        sam="mapped_to_assembly/{sample}.sam",
        fai="metaspades_output/{sample}/contigs.fasta.fai"
    output:
        "mapped_to_assembly/{sample}.bam"
    shell:
        "samtools view -bt {input.fai} {input.sam} > {output}"
rule sort_bam2:
    input:
        "mapped_to_assembly/{sample}.bam"
    output:
        "mapped_to_assembly/{sample}_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"
rule index_bam:
    input:
        "mapped_to_assembly/{sample}_sorted.bam"
    output:
        "mapped_to_assembly/{sample}_sorted.bam.bai"
    shell:
        "samtools index {input}"
rule generate_assembly_stats:
    input:
        "mapped_to_assembly/{sample}_sorted.bam"
    output:
        "assembly_stats/{sample}_idxstats.txt"
    shell:
        "samtools idxstats {input} > {output}"
rule generate_count_table:
    input:
        "assembly_stats/{sample}_idxstats.txt"
    output:
        "assembly_stats/{sample}_counts.txt"
    shell:
        "python get_count_table.py {input} > {output}"

rule extract_aligned_reads_for_mapped_reads:
    input:
        bam="mapped_reads/{sample}_host_mapped.bam"
    output:
        conserved_reads="extracted_reads_cyB/{sample}_conserved_region.bam"
    params:
        region="conserved_gene_coordinates.bed"  # BED file with coordinates of the conserved region
    shell:
        "samtools view -b {input.bam} -L {params.region} > {output.conserved_reads}"

#metaspades_take two for mapped host reads for extracting cyb gene:)
# Convert sorted BAM to FastQ
rule bam_to_fastq_for_cyb:
    input:
    	"mapped_reads/{sample}_conserved_region.bam"
    output:
    	R1="filtered_reads_cyb/{sample}_host_cyb_R1.fastq",
        R2="filtered_reads_cyb/{sample}_host_cyb_R2.fastq"
    shell:
    	"samtools fastq -1 {output.R1} -2 {output.R2} {input}"

rule run_metaspades_take_two:
    input:
        R1="filtered_reads_cyb/{sample}_host_cyb_R1.fastq",
        R2="filtered_reads_cyb/{sample}_host_cyb_R2.fastq"
    output:
        contigs="metaspades_output_cyb/{sample}/contigs.fasta",
        scaffolds="metaspades_output_cyb/{sample}/scaffolds.fasta"
    params:
        spades_path="spades.py",  # Update with your path to SPAdes
        threads=48,  # Adjust as needed
        memory=132,  # Adjust as needed, in GB
        outdir=lambda wc: "metaspades_output_cyb/" + wc.sample  # Dynamically generate output directory path
    shell:
        """
        {params.spades_path} --meta -1 {input.R1} -2 {input.R2} -k 55,75,95 \
        -t {params.threads} -m {params.memory} -o {params.outdir}
        """

rule centrifuge_to_krona:
    input:
        "centrifuge_reports/{sample}_centrifuge_report.tsv"
    output:
        "krona_charts/{sample}_krona.html"
    shell:
        "ktImportTaxonomy -q 1 -t 2 -m 7 {input} -o {output}"


rule index_bam2:
    input:
        "filtered_reads/{sample}_unmapped.bam"
    output:
        "filtered_reads_bam/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule filter_reads_by_taxon:
    input:
        report='centrifuge_reports/{sample}_centrifuge_report.tsv'
    output:
        ids='filtered_ids/{sample}_filtered_read_ids.txt'
    params:
        taxon_of_interest='248061'  # You can adjust this as needed
    script:
        "python filter_reads_by_taxon.py {input.report} {params.taxon_of_interest} {output.ids}"

# Rule to extract reads for each sample
rule extract_reads:
    input:
        ids='filtered_ids/{sample}_filtered_read_ids.txt',
	r1_fastq="filtered_reads/{sample}_nohost_R1.fastq",
        r2_fastq="filtered_reads/{sample}_nohost_R2.fastq"
    output:
        r1_fastq='extracted_reads/{sample}_R1_extracted_reads.fastq',
        r2_fastq='extracted_reads/{sample}_R2_extracted_reads.fastq'
    shell:
        """
        seqtk subseq {input.r1_fastq} {input.ids} > {output.r1_fastq}
        seqtk subseq {input.r2_fastq} {input.ids} > {output.r2_fastq}
        """

