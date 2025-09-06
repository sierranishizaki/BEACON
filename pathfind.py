import os
import sys
import time
import getopt

inputCache = 'cache'
outputCache = 'output'
try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=", "output="])
except getopt.GetoptError:
    print('python pathfind.py [-i INPUTCACHEDIR] [-o OUTPUTDIR]')
    sys.exit(2)
    
for opt, arg in opts:
    if opt in ('-h', '--help'):
        print("python pathfind.py [-i INPUTCACHE] [-o OUTPUTDIR]\n" +
              "INPUTCACHE defaults to 'cache', OUTPUT defaults to 'output'")
        sys.exit()
    elif opt in ("-i", "--input"):
        inputCache = arg
    elif opt in ("-o", "--output"):
        outputCache = arg

if not os.path.exists(outputCache):
    os.mkdir(outputCache)
if len(os.listdir(outputCache)) != 0:
    print(outputCache, "was not an empty directory, please delete files in", outputCache, "or use a different directory as a cache")

allProt = []

# this PQ a bit too fancy for the simple pathfinding algorithm that
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

tic = time.time()
counter = 0

bigCache = {}
with open(inputCache + "/allProteins.txt", 'r') as allPFile:
    for p in allPFile:
        if p[-1] == "\n":
            p = p[:-1]
        with open(inputCache + "/" + p, 'r') as f:
            bigCache[p] = []
            for d in f:
                bigCache[p].append(d.split()[0])

for source in sorted(bigCache.keys()):
    counter += 1
    if counter % 100 == 0:
        print(counter, "\t", time.time()-tic)

    distances = {source: 0}
    pq = [[source]]
    innerCounter = 0
    while sizePQ(pq) > 0:
        innerCounter += 1
        #if innerCounter % 1000 == 0:
        #    print(source, innerCounter, sizePQ(pq), time.time()-tic, sep="\t")
        (curProt, curSteps) = popPQ(pq)
        # if this protein doesn't have any interactions, skip it
        #   (this should only happen for proteins at step 0)
        #if not os.path.exists("cache/" + curProt):
        #    continue
        # if we've already seen this protein, skip it
        #   (the distance part shouldn't be necessary due to priority queue,
        #    but it is cheap to check just in case) TODO: remove?
        #if curProt != source and curProt in distances and curSteps >= distances[curProt]:
        #    print("Prevented revisiting", curProt)
        #    continue

        for nextProt in bigCache[curProt]:
            nextStep = 1
            # if we already know a faster way to get to this protein, skip it
            #   before the and prevents backtracking
            #   TODO: verify if the distance stuff adds value
            if nextProt in distances:
                continue
            addPQ(pq, nextProt, curSteps + 1)
            distances[nextProt] = curSteps + 1
        f.close()
    with open(outputCache + "/" + source, "w") as f:
        # sort the dictionary by value, then print
        sorted_items = sorted(distances.items(), key=lambda item: item[1])
        for prot, steps in sorted_items:
            print(prot, steps, file=f)
