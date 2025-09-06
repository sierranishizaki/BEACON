import gzip
import sys
import getopt
import os
import time
from pathlib import Path

filein = ''
threshold = 900
cache_dir = 'cache'

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:t:c:", ["help", "input=", "threshold=", "cache="])
except getopt.GetoptError:
    print('python initialParse.py -i INPUTFILE [-t INTEGER0-1000] [-c CACHEDIR]')
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        print("python initialParse.py -i INPUTFILE [-t INTEGER0-1000] [-c CACHEDIR]\n" +
              "INPUTFILE can be .txt or .txt.gz, threshold defaults to 900, cache defaults to 'cache'")
        sys.exit()
    elif opt in ("-i", "--input"):
        filein = arg
    elif opt in ("-t", "--threshold"):
        threshold = int(arg)
    elif opt in ("-c", "--cache"):
        cache_dir = arg

if not os.path.exists(filein):
    print("Input file '" + filein + "' was not found")
    sys.exit()
        
tic = time.time()

allProt = set()

counter = 0
modVal = 1000000

if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)
if len(os.listdir(cache_dir)) != 0:
    print(cache, "was not an empty directory, please delete files in", cache_dir, "or use a different directory as a cache")

if filein[-4:] == '.txt':
    f = open(filein, 'r')
elif filein[-7:] == '.txt.gz':
    f = gzip.open(filein, 'rt')
else:
    print("Unrecognized file format for '" + filein + "', please use .txt or .txt.gz")
    sys.exit()
    
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
            with open(cache_dir+"/"+prot1, "a") as w:
                w.write(prot2 + "\t1\n")
            with open(cache_dir+"/"+prot2, "a") as w:
                w.write(prot1 + "\t1\n")
        
with open(cache_dir + "/allProteins.txt", "w") as f:
    for p in sorted(allProt):
        print(p, file=f)
        Path(cache_dir+"/"+p).touch()

print(time.time()-tic)
