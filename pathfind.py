import os
import time

allProt = []


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
with open("allProteins.txt", 'r') as allPFile:
    for source in allPFile:
        counter += 1
        if counter % 100 == 0:
            print(counter, "\t", time.time()-tic)
        source = source[:-1] # trim off newline

        distances = {source: 0}
        pq = [[source]]
        while sizePQ(pq) > 0:
            (curProt, curSteps) = popPQ(pq)
            # if this protein doesn't have any interactions, skip it
            #   (this should only happen for proteins at step 0)
            if not os.path.exists("cache/" + curProt):
                continue
            # if we've already seen this protein, skip it
            #   (the distance part shouldn't be necessary due to priority queue,
            #    but it is cheap to check just in case) TODO: remove?
            #if curProt != source and curProt in distances and curSteps >= distances[curProt]:
            #    print("Prevented revisiting", curProt)
            #    continue

            # if we've already computed distances from this protein,
            #   use that work
            if os.path.exists("cache2/"+curProt):
                f = open("cache2/" + curProt, 'r')
            else:
                f = open("cache/" + curProt, 'r')

            for line in f:
                [nextProt, nextStep] = line.split()
                nextStep = int(nextStep)            
                # if we already know a faster way to get to this protein, skip it
                #   before the and prevents backtracking
                #   TODO: verify if the distance stuff adds value
                if nextProt in distances and distances[nextProt] < curSteps + nextStep:
                    continue
                addPQ(pq, nextProt, curSteps + nextStep)
                distances[nextProt] = curSteps + nextStep
            f.close()
        with open("cache2/" + source, "w") as f:
            # sort the dictionary by value, then print
            sorted_items = sorted(distances.items(), key=lambda item: item[1])
            for prot, steps in sorted_items:
                print(prot, steps, file=f)
