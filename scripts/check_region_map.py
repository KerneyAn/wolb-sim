import os

wmel_contig = "NC_002978.6"
wri_contig = "NC_012416.1"

d  = {"wmel":  wmel_contig, "wri": wri_contig}



# Function to parse filename and extract positions
def parse_filename(filename):
    parts = filename.split("-")
    if "-+-" in filename:
        donor_parts = parts[0:6]
        recipient_parts = parts[7:]
        
    else:
        donor_parts = parts[0:7]
        
        recipient_parts = parts[8:]

    donor_chrom = d[donor_parts[0]]
    recipient_chrom = d[recipient_parts[0]]
    donor_start, donor_end = map(int, donor_parts[4:6])
    recipient_start, recipient_end = map(int, recipient_parts[1:3])

    # donor_chrom = d[parts[0]]
    # recipient_chrom = d[parts[3]]
    # donor_start, donor_end = map(int, parts[4:6])
    # recipient_start, recipient_end = map(int, parts[8:10])

    return donor_chrom, donor_start, donor_end, recipient_chrom, recipient_start, recipient_end

# Function to read BED file and return mappable regions as a dict
def read_bed_file(bed_file):
    mappable_regions = {}
    with open(bed_file, 'r') as file:
        for line in file:
            chrom, start, end = line.strip().split()[:3]
            if chrom not in mappable_regions:
                mappable_regions[chrom] = []
            mappable_regions[chrom].append((int(start), int(end)))
    return mappable_regions

# Function to calculate mappable length in a given range
def calculate_mappable_length(chrom, start, end, mappable_regions):
    print(chrom)
    if chrom not in mappable_regions:
        return 0
    l = [(region_start, region_end) for region_start, region_end in mappable_regions[chrom] if region_start < end and region_end > start]
    # l = [(region_start, region_end) for region_start, region_end in mappable_regions[chrom] if region_start < end and region_start >= start and region_end > start and region_end <= end]
    print("start =", start,"end = ", end, "mappable regions:", l)  #TESTING

    num_map_bases = 0
    for sec in l:
        # print("region", sec) #TESTING
        num_map_bases += sec[1] - sec[0]
        # print("num_map_bases:", num_map_bases) #TESTING
    length = end - start
    perc_map_area = num_map_bases/ length


    print("length:", length) #TESTING
    print("num_map_bases:", num_map_bases) #TESTING
    print("percentage", perc_map_area) #TESTING
    print("\n") #TESTING

    # return(start, end, l, num_map_bases, length, perc_map_area )
    return sum(min(end, region_end) - max(start, region_start) for region_start, region_end in mappable_regions[chrom] if region_start < end and region_end > start)

# WHY!!!!!!!!!!!!!!!!!!!
# start(48398) end(48698) 
# total length: 300
# (47602, 49289) checking  this area for mappability but its larger than the total lengthhhhhhh, also its not even in the range
# this messes up my percentage of map area 5.623333333333333 but maybe this is okay because it just tells us that this mappable regin is larger than the the section we're looking at


# Main function to process recombinant genomes
def process_recombinant_genomes(bed_file, filenames):
    mappable_regions = read_bed_file(bed_file)
    results = []

    for filename in filenames:
        donor_chrom, donor_start, donor_end, recipient_chrom, recipient_start, recipient_end = parse_filename(filename)
        inserted_mappable_length = calculate_mappable_length(donor_chrom, donor_start, donor_end, mappable_regions)
        
        # Calculate surrounding mappable length including upstream and downstream 300bp regions
        upstream_start = recipient_start - 1000
        upstream_end = recipient_start
        downstream_start = recipient_end
        downstream_end = recipient_end + 1000

        print("up")
        upstream_mappable_length = calculate_mappable_length(recipient_chrom, upstream_start, upstream_end, mappable_regions)
        
        print("down")
        downstream_mappable_length = calculate_mappable_length(recipient_chrom, downstream_start, downstream_end, mappable_regions)
        
        surrounding_mappable_length = upstream_mappable_length + downstream_mappable_length

        results.append((filename, inserted_mappable_length, surrounding_mappable_length))
        results.append(upstream_mappable_length)
    print(results)
    return results

# Example usage
bed_file = '/home/kerney/wolb-recomb/simulations/testing_3_genome/results/mappability/0/min_1_map.bed'
file_names = [
    "wmel-into-wri-wmel-51488-74727-+-wri-48698-73288-+.fa",
    # "wmel-into-wri-wmel-462706-486423---wri-251266-278469-+.fa" 
    # Add other filenames as needed
]
# file_names = os.listdir(r"/home/kerney/wolb-recomb/data/simulated_genomes")
results = process_recombinant_genomes(bed_file, file_names)

# with open('file.txt', 'w') as f:
#     for result in results:
#         print(f"Filename: {result[0]}, Inserted Mappable Length: {result[1]}, Surrounding Mappable Length: {result[2]}")

for result in results:
     print(f"Filename: {result[0][0]}, Inserted Mappable Length: {result[0][1]}, Surrounding Mappable Length: {result[1]}")