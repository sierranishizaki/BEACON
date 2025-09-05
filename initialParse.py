import gzip
import os
import time
from pathlib import Path

tic = time.time()

allProt = set()

counter = 0
modVal = 1000000

if len(os.listdir("cache/")) != 0:
    print("cache was not empty - are you sure you didn't want to delete it first?")

with gzip.open("9606.protein.links.v12.0.txt.gz", 'rt') as f:
    for line in f:
        counter += 1
        if counter == 1:
            continue
        if counter%modVal == 0:
            print(counter//modVal, "\t", len(allProt), "\t", time.time()-tic)
        [prot1, prot2, score] = line.split()
        allProt.add(prot1)
        allProt.add(prot2)
        if prot1 < prot2 and int(score) >= 900:
            with open("cache/"+prot1, "a") as w:
                w.write(prot2 + "\t1\n")
        
with open("allProteins.txt", "w") as f:
    for p in sorted(allProt):
        print(p, file=f)
        Path("cache/"+p).touch()

print(time.time()-tic)
