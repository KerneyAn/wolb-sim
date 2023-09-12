from Bio import SeqIO
import random


# #wMel
# with open("/Users/kerney/RussellLab/RecombinationSims/wolb-sim/genomes/GCF_000008025.1_ASM802v1_genomic.fna") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         print(record.id)

# #wRi
# with open("/Users/kerney/RussellLab/RecombinationSims/wolb-sim/genomes/GCF_000022285.1_ASM2228v1_genomic.fna") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         print(record.id)

def faReader(faFilepath):
    tm = SeqIO.read(faFilepath, "fasta")
    return(str(tm.seq))


def Recombination(seqA, seqB, whichSeq, pos1, pos2):
    '''
    This is for site-specific recombination.
    Recombination(Sequence A, Sequence B, which Sequence(should be string i.e = "A"), Position1 of Sequence A, Position2 of Sequence B)
    '''
    seqList = [seqA, seqB]

    # Subtract 1 from all positions because python starts counting from 0
    recombPosone = tuple(map(lambda i, j: i - j, pos1, (1,0)))
    recombPostwo = tuple(map(lambda i, j: i - j, pos2, (1,0)))

    #Choosing the donor
    if whichSeq == "A":
        replaceSeq = seqA[recombPosone[0]:recombPosone[1]]
        recombSeq = seqB[:recombPostwo[0]] + replaceSeq + seqB[recombPostwo[1]:]
    else:
        replaceSeq = seqB[recombPostwo[0]:recombPostwo[1]]
        recombSeq = seqA[:recombPosone[0]] + replaceSeq + seqB[recombPosone[1]:]

    
    return(replaceSeq, recombPosone, recombPostwo, recombSeq)
    
# file path of sequences    
sequenceA = faReader("/Users/kerney/RussellLab/RecombinationSims/wolb-sim/genomes/GCF_000008025.1_ASM802v1_genomic.fna")
sequenceB = faReader("/Users/kerney/RussellLab/RecombinationSims/wolb-sim/genomes/GCF_000022285.1_ASM2228v1_genomic.fna")

# Sequence A (Pos1) and Sequence B (Pos2) homologous regions 
pos1 = (1,46036) 
pos2 = (1,46537) 

# #testing
# sequenceA = "GGTTTAAAGATCAAAATTAATGACAACAATTCTAAGGGTAAAGTTATGATACGGTATGATAATCCAAATGAATTGGATTTGATATTGAAAATTTTGAATAGGAAA"
# sequenceB = "TTAAACGCTTGACAAGCAAGATTAAGCTGCTTTTAATTGCAACTAACTTACGCTGCAAATGTTTAAGAAGTTTACTAAGCAGAAAAAAAGACAAAGAATCCCCGA"
# pos1 = (1,10) 
# pos2 = (1,12) 

# Choosing the donor
whichSequence = "A"


replaceSeq, recombPosone, recombPostwo, recombSeq = Recombination(sequenceA, sequenceB, whichSequence, pos1, pos2)


# print(recombSeq)

with open("test.fa", "w") as f:
    f.write(">test\n")
    f.write(recombSeq)