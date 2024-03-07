So i checked for the regions upstream and downstream by 1000 bases. I found that the mappable region upstream  found is larger than the length that I want to check.

So i made another l variable that will be stricter on the boundaries. But i dont think this is right because now it says there is no mappable region even though the  entire 1000 bases im checking is in the mappable region.

I found that the upstream of the inserted region it is much more mappable than the downstream. The single region is just larger. 

Here is a sample output

<!-- 
up
start = 47698 end =  48698 mappable regions: [(47602, 49289)]
region (47602, 49289)
num_map_bases: 1687
length: 1000
percentage 1.687


down
start = 73288 end =  74288 mappable regions: [(73573, 73764), (73796, 73985), (73990, 74288)]
region (73573, 73764)
num_map_bases: 191
region (73796, 73985)
num_map_bases: 380
region (73990, 74288)
num_map_bases: 678
length: 1000
percentage 0.678


Filename: wmel-into-wri-wmel-51488-74727-+-wri-48698-73288-+.fa, Inserted Mappable Length: 15860, Surrounding Mappable Length: 1678 
-->

I found the best map areas up down and at the insert. put them into bestMap.txt