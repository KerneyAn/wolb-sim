configfile: "config/config.yaml"


rule all:
    input:
        csv = "possible_recomb.txt"


rule merge_fa:
    input:
        re = config["recomb_genome_path"],
        wmel = config["wmel_genome_path"],
        wri = config["wri_genome_path"],
        dmel = config["dmel_genome_path"],
    output:
        out = "merged_genome/merged.fa"
    params:
        re = config["num_recomb"],
        wmel = config["num_wmel"],
        wri = config["num_wri"],
        dmel = config["num_dmel"],
    conda:
        "envs/env.yaml"
    script:
        "merge_fa.py"

rule sim_Genomes:
    input:
        merged = "merged_genome/merged.fa"
    output:
        out = ("simulated_reads/sim.bfast.fastq.gz"),
        read1 = "simulated_reads/sim.bwa.read1.fastq.gz",
        read2 = "simulated_reads/sim.bwa.read2.fastq.gz",
        mut = ("simulated_reads/sim.mutations.vcf"),
        mut_txt = ("simulated_reads/sim.mutations.txt")
    params:
        error_rate = config["sim_read_error_rate"],
        number_of_reads = config["num_reads"],
        out_name = "simulated_reads/sim"
    conda:
        "envs/env.yaml"
    shell:
        """
        dwgsim -N {params.number_of_reads} -1 150 -2 150 -H -y 0 -e {params.error_rate} -E {params.error_rate} -r 0 -F 0 {input.merged} {params.out_name}
        """

rule check_blocks:
    input:
        "simulated_reads/sim.bwa.read1.fastq.gz"
    output:
        txt = "possible_recomb.txt"
    params:
        ampRange = 150,
        recomb_genome = config["recomb_genome_path"]
    conda:
        "envs/env.yaml"
    script:
        "recom_check.py"


