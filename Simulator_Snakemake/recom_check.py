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
        ''' Parses the data from the headers into class variables'''
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

    
    def blockCheck(self, wmelDon, addRange=0):
        '''Checks for swapped blocks that are in range of homology'''

        count_in_range = 0
        count_out_of_range = 0
    
        if wmelDon == True:
            whichDonor = "wMel into wRi"
            # Check each element in the tuple and count how many are within or outside the specified range
            for value in (self.insertRange[0], self.insertRange[1]):
                if self.wriRange[0] <= value <= self.wriRange[1]:
                    count_in_range += 1
                else:
                    count_out_of_range += 1
        else:
            whichDonor = "wRi into wMel"
            for value in (self.insertRange[0], self.insertRange[1]):
                if self.wmelRange[0] <= value <= self.wmelRange[1]:
                    count_in_range += 1
                else:
                    count_out_of_range += 1
                
        # Check if one value is within the range and one value is outside the range
        if count_in_range == 1 and count_out_of_range == 1:
            return(True, whichDonor, self.head, (self.insertRange[0], self.insertRange[1]))
        else:
            return(False, whichDonor, self.head, (self.insertRange[0], self.insertRange[1]))
  

if __name__ == "__main__":
    # outfile = "checking.csv"
    outfile = Path(snakemake.output["csv"]) 
    addRan = snakemake.params["ampRange"]
    

    with open(outfile, 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        title = ["File Name", "Header", "Which Donor", "Range", "Total"]
        writer.writerow(title)

        for i in range(0,len(snakemake.input)):
            relName = snakemake.input[i]
            fullName = Path(__file__).parent / relName
        

            if "wmel-into-wri" in relName:
                wmelDonor = True
            elif "wri-into-wmel" in relName:
                wmelDonor = False
            
            thisPair = isPairRecomb()
            
            recombTotal = 0

            with fullName.open() as f:
                for line in f:
                    if line.startswith("@"):
                        header = line
                        thisPair.blockRanges(header)
                        inRange, whichDon, recomHead, insert = thisPair.blockCheck(wmelDonor, addRan)
                        if inRange == True:
                            info = [relName, recomHead, whichDon, insert]
                            writer.writerow(info)
                            recombTotal += 1


                spacer = [relName, "-", "-", "-", recombTotal]
                writer.writerow(spacer)           
                        

    
