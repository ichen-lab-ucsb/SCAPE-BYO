#Note: This is a python script to generate valid evolutionary pathways over the population of sequences
#present in a round of artificial selection. Priority is given to pathways with shorter steps and higher
#minimum sequence counts.

#Current input file requirements are "counts" files consisting of three lines of metadata 
#followed by one line per unique sequence in the pool in the following format: 
#sequences in the first column and counts (an integer number) in the second column. 

#Such files are produced by our Galaxy tools, currently available at the Chen Lab website:
# https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html

#by Abe Pressman
#contact: abe@engineering.ucsb.edu

#use as: python peak_pather.py counts_file output_file start_sequence end_sequence [optional arguments]

import Levenshtein
from heapq import heappop, heappush
import argparse

def main():
	parser = argparse.ArgumentParser(description="Find the shortest N evolutionary paths between two sequences in a population, using the A* algorithm")
	parser.add_argument("input", help="location/name of file containing sequence counts for the round to be iterated over")
	parser.add_argument("output", help="output file location/name")
	parser.add_argument("start_seq", help="path start sequence")
	parser.add_argument("end_seq", help="path end sequence")
	parser.add_argument("-i", "--in_type", default="counts", help="set input file type; default is 'counts', which assumes three lines of header data followed by lines of the format 'sequence count' where count is an integer")
	parser.add_argument("-n","--num_paths", type=int, default=1, help='number of best pathways to be generated with each run; default is 1')
	parser.add_argument("--max_length", type=int, default=0, help='maximum length of path searched before giving up; defaults to twice the length of the starting sequence')
	parser.add_argument("--min_count", type=int, default=2, help='minimum count of sequences searched (the program discards any sequences of lower count); default is 2 (setting lower than 2 may increase runtime dramatically)')
	parser.add_argument("--max_step", type=int, default=1, help='maximum step size allowed by search paths; default is 1')
	parser.add_argument("-d","--dist_type", default='edit', help='distance metric used to find shortest pathway; defaults to "edit" but "hamming" is also allowed')
	parser.add_argument("-p","--track_progress", action='store_true', help='if this flag is enabled, terminal will output updates on how much progress the code has made')

	args=parser.parse_args()

	maxPath = args.max_length
	if args.max_length == 0:
		maxPath = 2*len(args.start_seq) #set our maximum path length   

	counts = readSeqs(args.input, args.in_type, args.min_count, args.track_progress) #generate a dictionary of all sequences present and their count

	bestPaths = astar(args.start_seq, args.end_seq, counts, maxPath, args.dist_type, args.max_step, args.num_paths, args.track_progress) #find the N best pathways

	with open(args.output,'w') as fo:

		fo.write('These are the best ' +  str(args.num_paths) + ' paths from ' + str(args.start_seq) + ' to ' + str(args.end_seq) + '\n')
		fo.write('of maximum length ' + str(maxPath) + ' and maximum step size (' + str(args.dist_type) + ') ' + str(args.max_step) + '\n')
		fo.write('using sequences of minimum abundance ' + str(args.min_count) + '\n')
        
		for path in bestPaths:
			
			initDist = path[0] - path[6][0]
			
			fo.write('\nhere is a path\n')
			fo.write('step #,sequence,step size,total distance,sequence count\n')
			fo.write (str(path[3]) + ',' + str(path[4]) + ',' + str(path[1]) + ',' + str(initDist) + ',' + str(counts[path[4]]) + '\n')

			for step in path[6]:
				fo.write (str(step[3]) + ',' + str(step[4]) + ',' + str(step[1]) + ',' + str(step[0]) + ',' + str(step[2]) + '\n')

# Create a hashmap keyed to all sequences present in the population, returning a dictionary of all sequences we want to search over, and the number of times they appear in that round's sequencing data
def readSeqs(loc, fileType='counts', minCount=2, trackProgress=False):
    #fileType: 'counts' refers to a list of sequences followed by count numbers, with three lines of header information
    #min_count is the minimum count threshold accepted
    
    allSeqCounts = {}
    
    with open(loc) as f:
        if fileType == 'counts':
            line0 = next(f)
            next(f)
            next(f)            
            
            uniqueSplit = [elem for elem in line0.strip().split()]
            uniques = int(uniqueSplit[-1]) #total number of unique sequences present in this round
            
            i = 0
                        
            for lineRead in f:
                line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
           
                if trackProgress:
                    i += 1
                    if i % 100000 == 0:
                        print "Read " + str(i) + " of " + str(uniques) +" sequences from counts file"
                
                if int(line[1]) >= minCount:
                    allSeqCounts[line[0]] = int(line[1])
                          
    return allSeqCounts
    #a dictionary of all sequences we want to search over, and their abundance 

#A naive method to find all possible connections from each particular sequence to its neighbors. "Tricks" are added to significantly speed this step up
def greedyEdgesTricks(seq1, seqCounts, distanceType='edit', maxStep=1):
    
    dist = 0
    edgeList = []
    
    if distanceType == 'edit':
        
        len1 = len(seq1)
        len2 = 0
        adjDist = 0
        sublen = len1/2
        
        for seq2 in seqCounts:
            len2 = len(seq2)     
            adjDist = maxStep + abs(len1-len2)
                                 
            if Levenshtein.distance(seq1[:sublen], seq2[:sublen]) <= adjDist:
                
                #Not worth calculating the distance between the full strings if the distance between half-strings is
                #sufficient to fall over the threshhold of max_step + length difference. This speeds up the method
                #approximately four-fold
                
                dist = Levenshtein.distance(seq1, seq2)
                
                if dist <= maxStep:
                    
                    if dist > 0:
                        edgeList.append((seq2, dist))
                        #if seq2 is close enough to seq1 (and not identical), add it to the list of valid edges
             
    if distanceType == 'hamming':
        
        #similar method, for hamming distance
        
        for seq2 in seqCounts:
            dist = Levenshtein.hamming(seq1, seq2)
            if dist <= maxStep:
                if dist > 0:
                    edgeList.append((seq2, dist))
                
    return edgeList

#main function for generating pathways 
def astar(startSeq, endSeq, fullCounts, maxPath, distanceType='edit', maxStep=1, numPaths = 1, trackProgress=False):
    
    minPaths = [] #a heap of the shortest paths generated; will automatically keep paths sorted
    endDists = {} #a map of the shortest possible path distances from each sequence to the end sequence  
    
    if distanceType == 'edit':
               
        seqCounts = {} #a dictionary of counts for all sequences valid to iterate over 
        
        for seq in fullCounts:
                        
            if (Levenshtein.distance(seq, startSeq)+Levenshtein.distance(seq, endSeq)) <= maxPath:
                seqCounts[seq] = fullCounts[seq]
                
                endDists[seq] = Levenshtein.distance(seq, endSeq)
                
                #only valid to retain sequences whose distance to the start sequence + distance to the end point is overall lower than the maximum path length
                
    elif distanceType == 'hamming':
        for seq in seqCounts:
            endDists[seq] = Levenshtein.hamming(seq, endSeq)
  
    pathsFound = False #have we found a valid pathway yet?
    
    backPaths = {} #a map of the shortest back-path to reach each sequence
      
    #Each object stored in our "minPaths" heap has the format (traveled distance, maximum distance, minimum count, total steps, sequence, [back dist, this dist, this count, this step num, this sequence], list of paths backwards to the starting sequence)
    
    heappush(minPaths, (endDists[startSeq], 0, fullCounts[startSeq], 0, startSeq, [0, 0, fullCounts[startSeq], 0, startSeq], []))
    #we initialize the heap with a single tuple containing the starting sequence and an empty back-path list

    winners = [] #the list of successful paths we've found
    visited = {} #a dictionary of all the sequences we've visited so far, and how many times/along how many different pathways we've visited them
    
    edgeSet = {} #a dictionary of all the sequences whose neighbors we've calculated so far
    
    successLength = 0

    i = 0
    
    while(pathsFound == False): #keep iterating until we find all the paths we need
        
        thisPath = heappop(minPaths)
        
        thisPathDist = thisPath[0]
        
        thisMaxDist = thisPath[1] #maximum size of any step on the path being investigated
        
        thisMinCount = thisPath[2] #maximum sequence count of any step on the path being investigated
        
        thisTotalSteps = thisPath[3] #total steps on the path being investigated. More steps is better (for the same total length), because it means smaller average steps
        newTotalSteps = thisTotalSteps + 1
        
        thisSeq = thisPath[4] #what sequence have we reached so far?
        
        newNode = False #has this sequence been found before?

        if trackProgress:
            i += 1
            if i % 100 == 0:
                print str(i) + " edges traveled"
                print str(len(minPaths)) + " edges in queue"
                print "distance traveled: " + str(thisPathDist)
                print "shortest possible remaining distance: " + str(thisMaxDist)
                print " " 
        
        if thisSeq in visited:
            if visited[thisSeq] < numPaths: #not worth looking at the sequence "node" we've just reached if we've already reached it at least N times along shorter paths
                visited[thisSeq] += 1
                newNode = True
            else:
                edgeSet[thisSeq] = 0 #if it has been reached enough times, we nullify the set of edges to visit
        else:
            visited[thisSeq] = 0
            newNode = True
            edgeSet[thisSeq] = greedyEdgesTricks(thisSeq, seqCounts, distanceType, maxStep) #only worth finding neighboring sequences for the ones we actually visit
        
        newPath = [thisPath[5]] + thisPath[6] #what is the new pathway we travel along to reach this path (including it)
        
        thisBackDist = thisPath[5][0]
                    
        if newNode:
                            
            for edge in edgeSet[thisSeq]:
                #if we've already visited this next point on a shorter path, it's not an edge worth checking
                #otherwise, generate a list for each neighboring sequence "node" based on the "edge" it has to travel along
            
                newSeq = edge[0]
                newBackDist = edge[1] + thisBackDist
                                    
                newMinCount = min(thisMinCount,fullCounts[newSeq])
                newMaxDist = max(edge[1],thisMaxDist)
                
                heappush(minPaths, (endDists[newSeq]+newBackDist, newMaxDist, newMinCount, newTotalSteps, newSeq, [newBackDist, edge[1], fullCounts[newSeq], newTotalSteps, newSeq],newPath))
                #add a new tuple to the heap, of format (minimum possible heuristic distance, sequence, step distance, path backwards); the heap will automatically sort by shortest distance

        if thisSeq == endSeq:
            if len(winners)<numPaths:
                successLength = thisBackDist #once we reach the end sequence, we consider all possible paths of the same length
            else:
                if thisBackDist > successLength: #any paths longer than the shortest possible path are discarded
                    pathsFound = True
            winners.append(thisPath)
            
    rankedWinners = sorted(winners, key = lambda x: (x[0], x[1], -x[3], -x[2]))
    #sort best paths by shortest total path, then shortest maximum step size, then highest total steps (equivalent to shortest average step), then highest count; we only take the N best of these
    
    bestWinners = rankedWinners[0:numPaths]
    
    return bestWinners
   
if __name__ == "__main__":
    main()
