samples = ['H5P1']

rule all:
    input:
        expand("data/{sample}.fwd.fastq.gz", sample=samples),
        expand("data/{sample}.rev.fastq.gz", sample=samples),
        expand("data/{sample}.clumped.fq.gz", sample=samples),
        expand("data/{sample}.filtered_by_tile.fq.gz", sample=samples),
        expand("data/{sample}.trimmed.fq.gz",sample=samples),
        expand("data/{sample}.filtered.fq.gz", sample=samples),
        expand("data/{sample}.filtered.tadpole.fq.gz", sample=samples),
        expand("data/{sample}.final.rev.fq.gz", sample=samples),
        expand("data/{sample}.final.fwd.fq.gz", sample=samples),
        expand("data/{sample}.fastq.lrd.gz", sample=samples),
        expand("data/porechop.lrd.{sample}.fq", sample=samples),
        expand("data/longread_filtered.{sample}.fastq", sample=samples),
        expand("data/{sample}.aligned", sample=samples),
        expand("data/{sample}.good_processed.fasta", sample=samples),
        expand("results/nucleotides_prodigal/{sample}.fa", sample=samples),
        expand("results/crt_output/{sample}", sample=samples),
        expand("results/pilercr_out/{sample}", sample=samples),
        expand("results/prokka{sample}", sample=samples),
        expand("results/cmscan_out/{sample}", sample=samples),
        expand("results/rpsblast_output{sample}", sample=samples),
        expand("results/hmmer_pfam.out{sample}", sample=samples),
        expand("results/hmmer_TIGRfam.out{sample}", sample=samples),
        expand("results/tree{sample}", sample=samples)

### PRE-PROCESSING FOR SHORT READS ###

rule get_reads:                          ### When a URL is used to obtain files ###
    output:
        fwd="data/{sample}.fwd.fastq.gz",
        rev="data/{sample}.rev.fastq.gz"
     params:
        fwd_url=get_fwd_url,
        rev_url=get_rev_url
    shell:
        "wget -O {output.fwd} {params.fwd_url} ; wget -O {output.rev} {params.rev_url}"


rule run_bbmap_clumpify:
    input:
        raw_fwd="data/{sample}.fwd.fastq.gz",
        raw_rev="data/{sample}.rev.fastq.gz"
    output:
        "data/{sample}.clumped.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "clumpify.sh -eoom -Xmx200m -da in1={input.raw_fwd} in2={input.raw_rev} out={output}"

rule run_bbmap_filter_by_tile:
    input:
        rules.run_bbmap_clumpify.output
    output:
        "data/{sample}.filtered_by_tile.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "filterbytile.sh -eoom -Xmx200m -da in={input} out={output}"

rule run_bbmap_bbduk_remove_adapters:
    input:
        rules.run_bbmap_filter_by_tile.output
    output:
        "data/{sample}.trimmed.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "bbduk.sh -eoom -Xmx200m -da in={input} out={output} ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=50 ref=adapters ftm=5 ordered"

rule run_bbmap_bbduk_remove_artefacts:
    input:
        rules.run_bbmap_bbduk_remove_adapters.output
    output:
        "data/{sample}.filtered.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "bbduk.sh -eoom -Xmx200m -da in={input} out={output} k=31 ref=artifacts,phix ordered cardinality"

rule run_tadpole:
    input:
        rules.run_bbmap_bbduk_remove_artefacts.output
    output:
        "data/{sample}.filtered.tadpole.fq.gz"
    conda:
        "envs_qc.yaml"
    shell:
        "tadpole.sh -Xmx200m -eoom -da in={input} out={output} mode=correct k=62 ecc ecco merge=t"

rule split_reads:
    input:
        rules.run_tadpole.output
    output:
        fwd="data/{sample}.final.fwd.fq.gz",
        rev="data/{sample}.final.rev.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "reformat.sh -eoom -da in={input} out1={output.fwd} out2={output.rev}"


rule get_longreads:                          ### to use for URL to obtain reads ###
    output:
        lrd="data/{sample}.fastq.lrd.gz"
    params:
        long_url=get_long_url
    shell:
        "wget -O {output.lrd} {params.long_url}"

### REMOVING ADAPTERS FROM LONG READS ###
rule run_porechop:
    input:
        "data/{sample}.lrd.gz"
    output:
        "data/porechop.lrd.{sample}.fq.gz"
    conda:
        "envs/env_qc.yaml"
    shell:
        "porechop -i {input} --format fastq.gz -o {output}"

rule unzip:           ### must be unzipped so it can be input into prinseq ###
    input:
        rules.run_porechop.output
    output:
        "data/porechop.lrd.{sample}.fq"
    shell:
        "gunzip {input} > {output}"

rule run_longread_prinseq:
    input:
        rules.unzip.output
    output:
        "data/longread_filtered.{sample}.fastq"
    params:
        prinseq_file = "data/longread_filtered.{sample}"
    shell:
        "perl tools/prinseq-lite-0.20.4/prinseq-lite -fastq {input} -lc_method dust -lc_threshold 70 -out_good {params.prinseq_file}"

### HYBRID ASSMEBLY for the short and long reads ##
rule run_unicylcer:
    input:
        fwd = rules.split_reads.output.fwd,
        rev = rules.split_reads.output.rev,
        long = rules.run_longread_prinseq.output
    output:
        "data/{sample}.aligned"
    conda:
        "envs/unicycler_env.yaml"
    shell:
        "unicycler -1 {input.fwd} -2 {input.rev} --long {input.long} -o {output}"

rule run_prinseq:
    input:
        rules.run_unicylcer.output
    output:
        "data/{sample}.good_processed.fasta"
    conda:
        "envs/env.yaml"
    params:
        id_map = report("seq_id_map{sample}", category="ID Map"),
        good = report("{sample}.good_processed", category="Aligned and Trimmed Sequence")
    shell:
        "perl tools/prinseq-lite-0.20.4/prinseq-lite -fasta {input} -out_good {params.good} -out_format 1 -min_len 150 -noniupac -seq_id_mappings {params.id_map} -seq_id "

rule second:          ### Needed due to input for following rules is a child to another output ###
    input:
        rules.run_unicylcer.output
    output:
        "results/{sample}"
    shell:
        "cat {input}/assembly.fasta > {output}"

### FEATURE PREDICTION ###
rule run_prodigal:
    input:
      rules.second.output
    output:
      gbk_files = report("results/annotated_reads_prodigal/{sample}.gbk",  category="Feature Prediction"),
      protein_files = report("results/protein_translations_prodigal/{sample}.fa",  category="Feature Prediction"),
      nucleotide_files = report("results/nucleotides_prodigal/{sample}.fa",  category="Feature Prediction")
    conda:
      "envs/prodigal_env.yaml"
    shell:
      "prodigal -i {input} -o {output.gbk_files} -a {output.protein_files} -d {output.nucleotide_files} -m -c"

rule run_crt:
    input:
        rules.second.output
    output:
        report("results/crt_output{sample}", category="Feature Prediction")
    conda:
        "envs/crt_env.yaml"
    shell:
        "java -cp tools/CRT1.2-CLI.jar crt -minRL 20 -maxRL 50 -minSL 20 -maxSL 60 -searchWL 7 {input} {output}"

rule run_pilercr:
    input:
        rules.second.output
    output:
        report("results/pilercr_out{sample}", category="Feature Prediction")
    conda:
        "envs/pilercr_env.yaml"
    shell:
        "tools/pilercr1.06/pilercr -in {input} -out {output} -minarray 5  -mincons 0.75 -maxspacer 100 -minid 0.90"

rule run_prokka:
    input:
        rules.second.output
    output:
        report("results/prokka{sample}", category="Feature Prediction")
    conda:
        "envs/prokka_env.yaml"
    shell:
        "prokka {input} --outdir {output}"

rule run_cmscan:
    input:
        rules.second.output
    output:
        cmscan_tblout = report("results/cmscan_out{sample}", category="Feature Prediction")
    conda:
        "envs/cmscan_env.yaml"
    shell:
        "cmscan --tblout {output.cmscan_tblout} --fmt 2 databases/Rfam/Rfam.cm {input}"


### FUNCTIONAL ANNOTATION ###
rule run_rpsblast:
    input:
        rules.run_prodigal.output.protein_files
    output:
        report("results/rpsblast_output{sample}", category="Functional Annotation")
    conda:
        "envs/rpsblast_env.yaml"
    shell:
        "rpsblast -db databases/Cog/Cog -query {input} -out {output} -evalue 1e-2"

rule run_hmmer_pfam:
    input:
        rules.run_prodigal.output.protein_files
    output:
        report("results/hmmer_pfam.out{sample}", category="Functional Annotation")
    conda:
        "envs/env.yaml"
    shell:
        "hmmsearch -o {output} --cut_ga databases/Pfam/Pfam-A.hmm {input}"

rule run_hmmer_TIGRfam:
    input:
        rules.run_prodigal.output.protein_files
    output:
        report("results/hmmer_TIGRfam.out{sample}", category="Functional Annotation")
    conda:
        "envs/env.yaml"
    shell:
        "hmmsearch -o {output} --cut_nc databases/TIGRfam/TIGRFAMs_15.0_HMM.LIB {input}"

rule run_interproscan:
    input:
        rules.run_prodigal.output.protein_files
    output:
        "results/interproscan_output{sample}"
    shell:
        "./interproscan.sh -i {input} -o {output}"

rule run_signalp:
    input:
        rules.second.output
    conda:
        "envs/env.yaml"
    shell:
        "tools/./signalp {input}"

rule run_tmhmm:
    input:
        rules.second.output
    conda:
        "envs/env.yaml"
    shell:
        "tools/tmhmm-2.0c/tmhmm {input}"

### 16s rRNA phylogenetic tree ###

rule run_barrnap:
    input:
        "data/{sample}.fa"
    output:
        rna_fasta = "results/rna{sample}.fa",
        rna_gff = "results/rna{sample}.gff"
    conda:
        "envs/env.yaml"
    shell:
        "barrnap {input} {output.rna_gff} --outseq {output.rna_fasta}"

rule extract_16srrna:
    input:
        gff = rules.run_barrnap.output.rna_gff,
        fasta = "data/{sample}.fa"
    output:
        "data/{sample}.fa-16S.fna"
    conda:
        "envs/env.yaml"
    shell:
        "./16s.sh {input.gff} {input.fasta}"

rule make_database:
    input:
        "databases/current_Bacteria_aligned.fa"
    output:
        "databases/reference{sample}"
    conda:
        "envs/env.yaml"
    shell:
        "makeblastdb -in {input} -title reference -dbtype nucl -out {output}"

rule run_blast:
    input:
        ref = rules.make_database.output,
        query = rules.extract_16srrna.output
    output:
        "data/sequences.reference{sample}"
    conda:
        "envs/env.yaml"
    shell:
        "blastn -db {input.ref} -query {input.query} sequences.fasta -out{output}"

rule run_muscle:
    input:
        rules.run_blast.output
    output:
        "results/muscle_out{sample}"
    conda:
        "envs/env.yaml"
    shell:
        "muscle -in {input} -out {output}"

rule run_fasttree:
    input:
        rules.run_muscle.output
    output:
        "results/tree{sample}"
    conda:
        "envs/env.yaml"
    shell:
        "fasttree {input} > {output}"
