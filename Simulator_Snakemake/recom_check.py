# @wmel-51488-74727-+-wri-48698-73288-+_449490_449118_1_0_0_0_0:0:0_0:0:0_0/1

import pysam
import re
from pathlib import Path
import csv


class isPairRecomb():

    def __init__(self):
        self.wmelRange = ()
        self.wriRange = ()
        self.insertRange = ()
        self.head = ""

    def blockRanges(self, head):
        self.head = head

        patt1 = "\-[0-9]*\-[0-9]*-" # These are blocks that were swapped
        patt2 = '\+\_[0-9]*\_[0-9]*\_' # This the read pair of new insertion

        # First block is wmel/ second block is wri
        swapBlocks = re.findall(patt1, self.head)
        newInsert = re.findall(patt2, self.head)

        # indices 1-2 are the ranges
        wmelRan = swapBlocks[0].split("-")
        self.wmelRange = (int(wmelRan[1]), int(wmelRan[2]))
        wriRan = swapBlocks[1].split("-")
        self.wriRange = (int(wriRan[1]), int(wriRan[2]))

        insertRan = newInsert[0].split("_")
        self.insertRange  = (int(insertRan[1]), int(insertRan[2]))

    
    def blockCheck(self, wmelDon):
        count_in_range = 0
        count_out_of_range = 0

        # Check each element in the tuple and count how many are within or outside the specified range
        for value in (self.insertRange[0], self.insertRange[1]):
            if self.wriRange[0] <= value <= self.wriRange[1]:
                count_in_range += 1
            else:
                count_out_of_range += 1

        # Check if one value is within the range and one value is outside the range
        if count_in_range == 1 and count_out_of_range == 1:
            return((self.insertRange[0], self.insertRange[1]), "YES In range: " + self.head)
        else:
            return((self.insertRange[0], self.insertRange[1]), "NOT In range: " + self.head)
  

if __name__ == "__main__":
    # fileName = parseIn(snakemake.input[0])
    relName = 'simulated_reads/wmel-into-wri-wmel-51488-74727-+-wri-48698-73288-+_sim.bwa.read1.fastq'
    fullName = Path(__file__).parent / relName

    if "wmel-into-wri" in relName:
        wmelDonor = True
    elif "wri-into-wmel" in relName:
        wmelDonor = False
    else:
        pass
    
    thisPair = isPairRecomb()
    header = ""
    
    with open('profiles1.csv', 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        title = ["File Name", "Header", "Range" ]
        writer.writerow(title)
        with fullName.open() as f:
            for line in f:
                if line.startswith("@"):
                    header = line
                    thisPair.blockRanges(header)
                    insert, recomHead = thisPair.blockCheck(wmelDonor)
                    if "YES" in recomHead:
                        info = [relName, recomHead, insert]
                        writer.writerow(info)

        spacer = ["", "", ""]
        writer.writerow(spacer)           
                        

    
