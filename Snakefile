configfile: "config/config.yaml"

include: "rules/common.smk"

ERROR_VAR = 0

def find_files(wc):
    f = checkpoints.run_simpy.get().output[0]
    file_name = []
    with open("file_names.txt", "r") as file:
        file_name = file.read().splitlines()
    return(file_name)
 

# rule all:
#     input: 
#          lambda wc: expand("simulated_reads/{swapped_genome}_sim.{rn}", swapped_genome=find_files(wc),
#         rn=["bfast.fastq.gz","bwa.read1.fastq.gz", "bwa.read2.fastq.gz", "mutations.vcf", "mutations.txt"]) 

rule all:
    input:
        csv = "possible_recomb.csv"


checkpoint run_simpy:
    output:
        "file_names.txt"
    conda:
        "envs/env.yaml"
    script:
        "sim.py"

rule merge_fa:
    input:
        # re = lambda wc: expand("mixed_blocks/{swapped_genome}.fa", swapped_genome=find_files(wc))
        re = "mixed_blocks/{swapped_genome}.fa"
    output:
        out = "merged_reads/{swapped_genome}_merged.fa"
    params:
        howMany = [1, 1, 1, 1, 1, 1, 1, 1], #[mel, ri, host, re]
        groups = config["coverage_groups"],
        path = config["genome_path"]
    conda:
        "envs/env.yaml"
    script:
        "merge_fa.py"

rule sim_Genomes:
    input:
        merged = "merged_reads/{swapped_genome}_merged.fa"
    output:
        out = ("simulated_reads/{swapped_genome}_sim.bfast.fastq.gz"),
        read1 = "simulated_reads/{swapped_genome}_sim.bwa.read1.fastq.gz",
        read2 = "simulated_reads/{swapped_genome}_sim.bwa.read2.fastq.gz",
        mut = ("simulated_reads/{swapped_genome}_sim.mutations.vcf"),
        mut_txt = ("simulated_reads/{swapped_genome}_sim.mutations.txt")
    params:
        error_rate = ERROR_VAR,
        out_name = "simulated_reads/{swapped_genome}_sim"
    conda:
        "envs/env.yaml"
    shell:
        """
        dwgsim -N 200000000 -1 150 -2 150 -H -y 0 -e {params.error_rate} -E {params.error_rate} -r 0 -F 0 {input.merged} {params.out_name}
        """


# rule sim_Genomes:
#     input:
#         "mixed_blocks/{swapped_genome}.fa"
#     output:
#         out = temp("simulated_reads/{swapped_genome}_sim.bfast.fastq.gz"),
#         read1 = "simulated_reads/{swapped_genome}_sim.bwa.read1.fastq.gz",
#         read2 = "simulated_reads/{swapped_genome}_sim.bwa.read2.fastq.gz",
#         mut = temp("simulated_reads/{swapped_genome}_sim.mutations.vcf"),
#         mut_txt = temp("simulated_reads/{swapped_genome}_sim.mutations.txt")
#     params:
#         error_rate = ERROR_VAR,
#         out_name = "simulated_reads/{swapped_genome}_sim"
#     conda:
#         "envs/env.yaml"
#     shell:
#         """
#         dwgsim -N 10000 -1 150 -2 150 -H -y 0 -e {params.error_rate} -E {params.error_rate} -r 0 -F 0 {input} {params.out_name}
#         """

rule extractRead1:
    input:
        "simulated_reads/{swapped_genome}_sim.bwa.read1.fastq.gz"
    output:
        "simulated_reads/{swapped_genome}_sim.bwa.read1.fastq"
    shell:
        """
        gunzip {input}
        """


rule check_blocks:
    input:
        lambda wc: expand("simulated_reads/{swapped_genome}_sim.bwa.read1.fastq", swapped_genome=find_files(wc))
    output:
        csv = "possible_recomb.csv"
    params:
        ampRange = 150
    conda:
        "envs/env.yaml"
    script:
        "recom_check.py"

