# - we want simulate from "realistic" scenario of mixed genomes
# 	- can make fasta containing different numbers of genomes to replicate real infections
# 	- so for example, make a fasta with a 10-30 copies of wmel reference, 1 copy of wri, some copies/different recombs and host
# 	- simulate from this composite fasta
# 	- make a script that can make this composite with arguments to control how many recombs, how many wmels, how many wris, etc
# fasta file:
# ```
# >genome1 wmel
# actg
# >genome2 wmel
# actg
# >host
# actg
# > wri
# > recomb
# > recomb
# ```

def read_genome(file_path):
    """Reads genome sequences from a FASTA file and returns a dictionary of header:sequence."""
    sequences = {}
    with open(file_path, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    # Save the previous sequence
                    sequences[header] = sequence
                # Start a new sequence
                header = line[1:]  # Remove '>'
                sequence = ''
            else:
                sequence += line
        if header:
            # Save the last sequence
            sequences[header] = sequence
    return sequences

def write_merged_genome(output_path, genomes, counts):
    """Writes the merged genome sequences to a FASTA file with unique headers."""
    with open(output_path, 'w') as file:
        for genome_name, genome_dict in genomes.items():
            for header, seq in genome_dict.items():
                for i in range(1, counts[genome_name] + 1):
                    # Write unique header
                    file.write(f'>{genome_name}_{header}_{i}\n')
                    # Write genome sequence
                    file.write(f'{seq}\n\n')



def main():

    # Load the genome sequences
    genomes = {
        'recomb': read_genome(snakemake.input.re),
        'wmel': read_genome(snakemake.input.wmel),
        'wri': read_genome(snakemake.input.wri),
        'dmel': read_genome(snakemake.input.dmel)
    }

    # Get the counts for each genome
    counts = {
        'recomb': int(snakemake.params.re),
        'wmel': int(snakemake.params.wmel),
        'wri': int(snakemake.params.wri),
        'dmel': int(snakemake.params.dmel)
    }

    # Write the merged genome
    write_merged_genome(snakemake.output[0], genomes, counts)

if __name__ == "__main__":
    main()