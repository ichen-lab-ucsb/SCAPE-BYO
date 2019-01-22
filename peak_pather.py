#Note: the current script is a placeholder/work in progress, meant to be tinkered with
#during use/analysis. It cannot currently be run from the command line, and 
#contains a number of diagnostic outputs to the console window, as well as
#additional functionality not necessary for main use.
#(An example of use is included, commented out, as an example of use)

#A new version, cleaned up and with better commenting, is being worked on.
#This version is included for use in the meanwhile, as necessary.

import Levenshtein
from scipy.spatial import distance
from heapq import heappop, heappush
import time
from matplotlib.pyplot import step

def readPeakCenters(loc, fileType='landscapeID'):
    
    #fileType: 'landscapeID' refers to the identified-landscape files used originally when writing this code
    #We only care about the peak centers, not whole peaks
    
    peakCtrList = []
    
    with open(loc) as f:
        if fileType == 'landscapeID':
            
            new = False
            
            for lineRead in f:
                line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]

                if line != []:
                    line[len(line)-1] = line[len(line)-1].split('\n')[0]
                    #remove newline character from the end of each line
                    
                    if line[1] == 'peak':
                        new = True        
                    else:
                        if new:
                            new = False
                            peakCtrList.append(line[1])
                            
    return peakCtrList
    #a list of all peak centers we might want to search between



def readSeqs(loc, fileType='text', cutOff=3):
    #fileType: 'counts' refers to a list of sequences followed by count numbers, with three lines of header information
    
    allSeqCounts = {}
    
    with open(loc) as f:
        if fileType == 'text':
            
            i = 0
            
            for lineRead in f:
                line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
                
                i += 1
                
                if i%100000 == 0:
                    print i
                
                if int(line[1]) >= cutOff:
                    allSeqCounts[line[0]] = int(line[1])
                    
        
                    
    return allSeqCounts
    #a dictionary of all sequences we want to search over, and their abundance 



def greedyEdges(seq1, seqCounts, distanceType='edit', maxDist=1):
    #a naive way to find all possible connections from a particular sequence to its neighbors
    
    dist = 0
    edgeList = []

    print 0
    
    print time.time()

    testval1 = 0
    testval2 = 0
    
    if distanceType == 'edit':
        
        for seq2 in seqCounts:
            dist = Levenshtein.distance(seq1, seq2)
            
            testval1 += 1
            
            if dist <= maxDist:
                
                testval2 += 1
                
                if dist > 0:
                    edgeList.append((seq2, dist))
                
    if distanceType == 'hamming':
        
        for seq2 in seqCounts:
            dist = Levenshtein.hamming(seq1, seq2)
            if dist <= maxDist:
                if dist > 0:
                    edgeList.append((seq2, dist))
                
                    
    return edgeList

def greedyEdgesTricks(seq1, seqCounts, distanceType='edit', maxDist=1):
    #a naive way to find all possible connections from a particular sequence to its neighbors, using weird tricks to speed up
    
    dist = 0
    edgeList = []

    

    
    if distanceType == 'edit':
        
        len1 = len(seq1)
        len2 = 0
        adjDist = 0
        sublen = len1/2
        
        for seq2 in seqCounts:
            len2 = len(seq2)            
            adjDist = maxDist + abs(len1-len2)
                        
                
            if Levenshtein.distance(seq1[:sublen], seq2[:sublen]) <= adjDist:
                dist = Levenshtein.distance(seq1, seq2)
                
                
                if dist <= maxDist:
                    
                    if dist > 0:
                        edgeList.append((seq2, dist))
                        

                
    if distanceType == 'hamming':
        
        for seq2 in seqCounts:
            dist = Levenshtein.hamming(seq1, seq2)
            if dist <= maxDist:
                if dist > 0:
                    edgeList.append((seq2, dist))
                
                
    
    return edgeList


def smartEdges(seq1, seqCounts, distanceType='edit', maxDist=1):
#a way to grab all possible neighbors and see if they exist. Runs faster when maxDist is low.
        
    searchDict = {}
    searchDict[seq1] = 0
    edgeList = []
    tempSeq = ''
    
    for i in range(maxDist):
        for seq in searchDict:
            for pos in range(len(seq)):
                for nuc in ['A','C','T','G']:
                    tempSeq = seq[:pos] + nuc + seq[pos+1:]
                    
                    if tempSeq not in searchDict:
                        searchDict[tempSeq] = i+1
                        
                if distanceType == 'edit':
                    for nuc in ['A','C','T','G']:
                        tempSeq = seq[:pos] + nuc + seq[pos:]
                    
                        if tempSeq not in searchDict:
                            searchDict[tempSeq] = i+1

                    
                    tempSeq = seq[:pos] + seq[pos+1:]
                    if tempSeq not in searchDict:
                        searchDict[tempSeq] = i+1
            if distanceType == 'edit':
                for nuc in ['A','C','T','G']:
                    tempSeq = seq + nuc
                    if tempSeq not in searchDict:
                        searchDict[tempSeq] = i+1
                        
    for seq, dist in searchDict:
        if seq in seqCounts:
            edgeList.append((seq, dist))
            
    return edgeList


def astar(startSeq, endSeq, fullCounts, maxDetours, distanceType='edit', maxDist=1, multiPaths='smallestSteps', numPaths = 1, distTolerance = 0):
    
    minPaths = []
    endDists = {}
    
    test1 = 0
    test2 = 0
    
    if distanceType == 'edit':
        maxPath = Levenshtein.distance(startSeq, endSeq) + maxDetours
        
        seqCounts = {}
        
        for seq in fullCounts:
            test1 += 1
            
            if (Levenshtein.distance(seq, startSeq)+Levenshtein.distance(seq, endSeq)) <= maxPath:
                seqCounts[seq] = fullCounts[seq]
                test2 += 1
                
                endDists[seq] = Levenshtein.distance(seq, endSeq)
    elif distanceType == 'hamming':
        for seq in seqCounts:
            endDists[seq] = Levenshtein.hamming(seq, endSeq)
    
    print test1
    print test2
    
    pathsFound = False
    
    
    backPaths = {}
    
    

    
    if numPaths == 1:     
        
        #backPaths[startSeq] = [0, 0, 0, seqCounts[startSeq], '']
        #traveled distance, step distance, total steps, minimum count, path backwards
        
        heappush(minPaths, (endDists[startSeq], startSeq, 0, ''))
        #minimum possible heuristic distance, sequence, step distance, path backwards

        timepoint = 0
           
        while(pathsFound == False):
            
            print time.clock() - timepoint
            timepoint = time.clock()
            
            thisPath = heappop(minPaths)
            vertex = thisPath[1]
            previous = thisPath[3]
            
            thisDist = thisPath[2]
            if previous in backPaths:
                thisDist += backPaths[previous][0]
                
            thisMinimum = seqCounts[vertex]
            if previous in backPaths:
                thisMinimum = min(seqCounts[vertex], backPaths[previous][3])
                                  
            thisStepNo = 1
            if previous in backPaths:
                thisStepNo += backPaths[previous][2]
    
            print 'part1 ' + str(time.clock() - timepoint) + '\n'
            
            if vertex in backPaths:
            #if we've already gotten here once...
                
                comparePath = backPaths[vertex]
                
                if comparePath[0] == thisDist:
                    #two equidistant paths to the same point are an option (our new path CANNOT be shorter because of priority)
                
                    #compare paths to get to this point
                    if multiPaths == 'highestCount':
                    #choose the path with the highest minimum sequence count 
                        if thisMinimum > comparePath[3]:
                            backPaths[vertex] = [thisDist, thisPath[2], thisStepNo, thisMinimum, previous]
                    elif multiPaths == 'smallestSteps':
                        if thisStepNo > comparePath[2]:
                            backPaths[vertex] = [thisDist, thisPath[2], thisStepNo, thisMinimum, previous]
                    elif multiPaths == 'fewestSteps':
                        if thisStepNo < comparePath[2]:
                            backPaths[vertex] = [thisDist, thisPath[2], thisStepNo, thisMinimum, previous]
                    #for other options: decide how we want to handle multiple paths; same structure? different one?
            
                print 'part2 ' + str(time.clock() - timepoint) + '\n'

            else:
                
                backPaths[vertex] = [thisDist, thisPath[2], thisStepNo, thisMinimum, previous]
                #traveled distance, step distance, total steps, minimum count, path backwards
    
                
                testEdges = greedyEdgesTricks(vertex, seqCounts, distanceType, maxDist)
                
                for edge in testEdges:
                    if edge[0] not in backPaths:
                        #if we've already visited this next point on a shorter path, it's not an edge worth checking
                    
                        nextVertex = edge[0]
                        distToNext = edge[1] + thisDist
                        
                        heappush(minPaths, (endDists[nextVertex]+distToNext, nextVertex, edge[1], vertex))
                        #minimum possible heuristic distance, sequence, step distance, path backwards
    
                print 'part3 ' + str(time.clock() - timepoint) + '\n'

            
            print thisPath[0]-thisDist
            print thisPath[0]
            
            if endSeq in backPaths:
                if backPaths[endSeq][0] < thisPath[0]:
                    pathsFound += 1
        
            print '\ndoop'
            
            
        nextVertex = endSeq
        
        print 'Sequence\tTotalDist\tStepDist\tSeqCount'
        
        while nextVertex:
            vtx = backPaths[nextVertex]
            print nextVertex + '\t' + str(vtx[0]) + '\t' + str(vtx[1]) + '\t' + str(seqCounts[nextVertex]) + '\n'
            
            nextVertex = vtx[4]
            
    else:
        
        #backPaths[startSeq] = [0, 0, 0, seqCounts[startSeq], [], []]
        #(traveled distance, maximum distance, minimum count, total steps, sequence, [back dist, this dist, this count, this step num, this sequence], list backwards)
        
        heappush(minPaths, (endDists[startSeq], 0, fullCounts[startSeq], 0, startSeq, [0, 0, fullCounts[startSeq], 0, startSeq], []))

        winners = []
        visited = {}
        
        edgeSet = {}
        
        successLength = 0

        timepoint = 0

        i = 0

        while(pathsFound == False):
            
            i += 1
            
 
        
            
            thisPath = heappop(minPaths)
            
            thisPathDist = thisPath[0]
            
            thisMaxDist = thisPath[1]
            
            thisMinCount = thisPath[2]
            
            thisTotalSteps = thisPath[3]
            newTotalSteps = thisTotalSteps + 1
            
            thisSeq = thisPath[4]
            
            newNode = False
            
            
            if thisSeq in visited:
                if visited[thisSeq] < numPaths:
                    visited[thisSeq] += 1
                    newNode = True
                else:
                    edgeSet[thisSeq] = 0
            else:
                visited[thisSeq] = 0
                newNode = True
                edgeSet[thisSeq] = greedyEdgesTricks(thisSeq, seqCounts, distanceType, maxDist)

            
            newPath = [thisPath[5]] + thisPath[6]
            
            thisBackDist = thisPath[5][0]
                        

            
            if newNode:
                                
                for edge in edgeSet[thisSeq]:
                    #if we've already visited this next point on a shorter path, it's not an edge worth checking
                
                    newSeq = edge[0]
                    newBackDist = edge[1] + thisBackDist
                                        
                    newMinCount = min(thisMinCount,fullCounts[newSeq])
                    newMaxDist = max(edge[1],thisMaxDist)
                    
                    heappush(minPaths, (endDists[newSeq]+newBackDist, newMaxDist, newMinCount, newTotalSteps, newSeq, [newBackDist, edge[1], fullCounts[newSeq], newTotalSteps, newSeq],newPath))
                    #minimum possible heuristic distance, sequence, step distance, path backwards


            if i % 100 == 0:
                print 'time  ' + str(time.clock()-timepoint)
            
                print '\n\n'
                timepoint = time.clock()   
    
                print thisPath[0]-thisBackDist
                print thisPath[0]
                print len(winners)



            if thisSeq == endSeq:
                if len(winners)<numPaths:
                    successLength = thisBackDist
                else:
                    if thisBackDist > successLength:
                        pathsFound = True
                winners.append(thisPath)
                
        rankedWinners = sorted(winners, key = lambda x: (x[0], x[1], -x[3], -x[2]))
        
        bestWinners = rankedWinners[0:numPaths]
        
        for path in bestWinners:
            print 'here is a path\n'
            print path[5]
            print '\n'
            for step in path[6]:
                print step
                print '\n'
            
        
#Example case:

#counts = readSeqs('R5c-counts-trim.txt', cutOff=2)
#astar('CCACACTTCAAGCAATCGGTC', 'CTCTTCAATAATCGGTTGCGT', counts, 35, maxDist=3, numPaths=5)

    