import gzip
import sys
import getopt
import os
import time
from pathlib import Path

# these define the default values
filein = ''
threshold = 400
output_dir = 'output'
silent = False

tic = time.time()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hsi:t:o:",
                               ["help", "silent", "input=", "threshold=", "output="])
except getopt.GetoptError:
    print('python beaconPathfind.py -i INPUTFILE [-s] [-t NUM0-1] [-o OUTPUTDIR]')
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        print("python beaconPathfind.py -i INPUTFILE [-s/--silent] [-t NUM0-1] [-o OUTPUTDIR]\n" +
              "INPUTFILE can be .txt or .txt.gz, output defaults to 'output'\n" +
              "threshold can be 0-1 (3 significant figures) or an integer 0-1000, defaults to 0.4")
        sys.exit()
    elif opt in ("-i", "--input"):
        filein = arg
    elif opt in ("-s", "--silent"):
        silent = True
    elif opt in ("-t", "--threshold"):
        # threshold could either be in decimal form (how the papers talk about it)
        if float(arg) < 1:
            threshold = int(1000*float(arg))
        # or integer form (how the database encodes it, therefore how this code handles it)
        else:
            threshold = int(arg)
        if threshold > 1000:
            print("Threshold is too high to be reasonable!")
    elif opt in ("-o", "--output"):
        output_dir = arg

if not silent:
    print("Running with parameters:\ninput: '" + filein + "'\nthreshold: '" +
          str(threshold/1000) + "'\noutput: '" + output_dir + "'")
        
if not os.path.exists(filein):
    print("Input file '" + filein + "' was not found")
    sys.exit()

counter = 0
modVal = 1000000

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if len(os.listdir(output_dir)) != 0:
    print(output_dir, "was not an empty directory, please delete files in", output_dir, "or use a different directory for output")
    sys.exit()

if filein[-4:] == '.txt':
    f = open(filein, 'r')
elif filein[-7:] == '.txt.gz':
    f = gzip.open(filein, 'rt')
else:
    print("Unrecognized file format for '" + filein + "', please use .txt or .txt.gz")
    sys.exit()

if not silent:
    print("Reading protein interactions from", filein, "with threshold", threshold/1000, "columns are:\nRows processed (*" + str(modVal) + "), Unique Proteins, Time elapsed")

# this reads in a file in the format of:
# GENE1 GENE2 Integer(0-1000)
# and saves any lines where the Integer is >= Threshold
# as a direct interaction. This will then be used in the next
# step for pathfinding
bigCache = {}
for line in f:
    counter += 1
    if counter == 1:
        continue
    if counter%modVal == 0 and not silent:
        print(counter//modVal, "\t", len(bigCache), "\t", time.time()-tic)
    [prot1, prot2, score] = line.split()
    if prot1 not in bigCache:
        bigCache[prot1] = []
    if prot2 not in bigCache:
        bigCache[prot2] = []
    if prot1 < prot2 and int(score) >= threshold:
        bigCache[prot1].append(prot2)
        bigCache[prot2].append(prot1)

if not silent:
    print("Done processing input and building cache, took", time.time()-tic)

# this PQ is a bit too fancy for the simple pathfinding algorithm that
#   ended up being used for this, it could have just been a flat queue

# pq is a priority queue organized to be a list of lists
#   the priority is an integer (number of steps from source)
#   and determines which internal list the protein is stored in
def sizePQ(pq):
    s = 0
    for arr in pq:
        s += len(arr)
    return s

def addPQ(pq, prot, steps):
    while len(pq) <= steps:
        pq.append([])
    pq[steps].append(prot)

# require that sizePQ(pq) > 0
def popPQ(pq):
    for i in range(len(pq)):
        if len(pq[i]) > 0:
            return (pq[i].pop(), i)

counter = 0

if not silent:
    print("Building networks of proteins, writing to '" + output_dir + "/', columns are:\n" +
          "# proteins processed (of "+str(len(bigCache))+") Time Elapsed from start of whole script")

# for each protein, build a network from it
for source in sorted(bigCache.keys()):
    counter += 1
    if counter % 100 == 0 and not silent:
        print(counter, "\t", time.time()-tic)

    distances = {source: 0}
    pq = [[source]]
    while sizePQ(pq) > 0:
        (curProt, curSteps) = popPQ(pq)
        for nextProt in bigCache[curProt]:
            nextStep = 1
            # if we already know a faster way to get to this protein, skip it
            if nextProt in distances:
                continue
            addPQ(pq, nextProt, curSteps + 1)
            distances[nextProt] = curSteps + 1
    with open(output_dir + "/" + source, "w") as f:
        # sort the dictionary by value, then print
        sorted_items = sorted(distances.items(), key=lambda item: item[1])
        for prot, steps in sorted_items:
            print(prot, steps, file=f)
if not silent:
    print("Total time elapsed:", time.time()-tic)
