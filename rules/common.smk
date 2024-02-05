from pathlib import Path

# Get a path object for the directory to list
directory = config["genome_path"]

genome_directory = Path(directory).glob('*.fa')
genome_name = []

for file in genome_directory:
    genome_name.append(Path(file).stem)





