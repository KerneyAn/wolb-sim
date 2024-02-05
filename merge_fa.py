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
class fastaMerger():
    def __init__(self):
        x = ""
    def merger(self, inputFiles, countFiles, outFile):
        with open(outFile, 'a') as out_fasta:
            file_counter = 1
            for c, numFiles in enumerate(countFiles):
                for num in range(numFiles):
                    with open(inputFiles[c], 'r') as in_fasta:
                        content = in_fasta.read()
                        # Add a counter to the header names to avoid duplicates
                        lines = content.split('\n')
                        for i,line in enumerate(lines):
                            if line.startswith('>'): # Check if it's a header line
                                line = f'{line[:-1]}_{file_counter}'
                                file_counter += 1
                            out_fasta.write(line)
                        
                            if i < len(lines) and line.strip():  # Add newline unless it's the last line
                                out_fasta.write('\n')
                               
                    
                # Reset the file counter for the next input file
                file_counter = 1

def main():

    thisMerger = fastaMerger()

    path = snakemake.params["path"]
    howMany = snakemake.params["howMany"]  
    swapped = snakemake.input["re"]
    groups = snakemake.params["groups"] 


    inFiles = []
    for k, v in groups.items():
        for i in v:
            file = path + i + ".fa"
            inFiles.append(file)

    files = inFiles.copy()
   
    

    # for m in swapped:
    #     name = m.replace("mixed_blocks/", "merged_reads/")
    #     name = name.replace(".fa", "")
    #     out = name + "_merged.fa"
    #     files.extend([m])
    #     thisMerger.merger(files, howMany, out)
    #     files = inFiles.copy()

    name = swapped.replace("mixed_blocks/", "merged_reads/")
    out = name.replace(".fa", "_merged.fa")
    files.extend([swapped])
    thisMerger.merger(files, howMany, out)

    # # testing
    # mel =  "x_samp/mel.fa"
    # re =  "x_samp/re.fa"
    # ri =  "x_samp/ri.fa"
    # outfileName = "merged.fa"

    # ho1  = "x_samp/a.fa"
    # ho2  = "x_samp/b.fa"
    # ho3  = "x_samp/c.fa"
    # ho4  = "x_samp/d.fa"
    # ho5  = "x_samp/e.fa"

    # inFiles = [mel, ri, ho1, ho2, ho3, ho4, ho5, re]
    # howMany = [30, 1, 1, 1]
    # howMany.extend(howmany[4])
    # thisMerger = fastaMerger(outfileName)
    # thisMerger.merger(inFiles, howMany)

    # for file in infiles:
    #     file = Path(file)
    #     depths.append(parse(file, groups))




if __name__ == "__main__":
    main()

