#Note: This is a python script to generate k-seq data over all sequences present in a round of selection,
#using a set of "k-seq rounds" as additional data to calculate the catalytic activity of every sequence

#Current input file requirements are "counts" files consisting of three lines of metadata 
#followed by one line per unique sequence in the pool in the following format: 
#sequences in the first column and counts (an integer number) in the second column. 

#Such files are produced by our Galaxy tools, currently available at the Chen Lab website:
# https://labs.chem.ucsb.edu/chen/irene/Chen_lab_at_UCSB/Galaxy_Tools.html

#by Abe Pressman
#contact: abe@engineering.ucsb.edu

#use as: python kseq_tools_v01.py start_round kseq_rounds output_file normalization_list substrate_concs rounds_to_average rounds_to_error 

import numpy as np
from scipy.optimize import curve_fit
import Levenshtein
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate catalytic kinetics for a population of sequences, using the k-seq methodology. Takes a k-seq 'start round' and a list of additional rounds; gives predicted constants A and k*t for catalysis following [surviving fraction] = A(1-Exp(-k*[S]*t)). If generally confused over use, see example use case given on github")
    parser.add_argument("start_round", help="location/name of file containing sequence counts for the pre-k-seq population")
    parser.add_argument("kseq_rounds", help="location/name of file containing a list of additional filenames, each of which contains sequence counts for a post-k-seq population")
    parser.add_argument("output", help="output file location/name")
    parser.add_argument("normalization_list", help="list of normalization factors for each round, starting with the start round. One value per line. For the start round, should typically be set equal to 1/[total amount of DNA/RNA/protein] present at start of k-seq selection rounds, and for all other rounds should be 1/[amount of DNA/RNA/protein] remaining after selection step.")    
    parser.add_argument("substrate_concs", help="file containing list of substrate concentrations for each set of kseq rounds")
    parser.add_argument("rounds_to_average", help="file containing comma-separated lists of of sets of rounds to average together for fits  (e.g. 1,2,3/4,5,6 will average the abundances of rounds 1,2, and 3, and then also average 4,5, and 6). If only a single replicate was carried out, then 1/2/3 etc. If you just want abundance without kseq calculations, leave blank. See example use case if confused")    
    parser.add_argument("rounds_to_error", help="file containing comma-separated lists of of sets of rounds used for each replicate. If no replicates, leave blank (and st.dev will be based on goodness of fit instead). See example use case if confused")    
    parser.add_argument("-s","--search_set", nargs='+', default='all', help="the set of sequences to search over. By default, set to 'all,' generating k-seq data for all sequences in the start round. If set to 'center [center_sequence] [distance]' (requires three arguments, the second being a sequence and the third being an integer), will only generate kseq data over all sequences within a fixed distance of the center. If set to 'list [sequence_list_location], will only generate data over the sequences listed, one on each row, in a file at the given location")
    parser.add_argument("-v","--verbose", action='store_true', help='if this flag is enabled, output file will include data on sequence concentration at every kseq round, resulting in a larger output file')
    parser.add_argument("-i", "--in_type", default="counts", help="set input file type; default is 'counts', which assumes three lines of header data followed by lines of the format 'sequence count' where count is an integer")
    parser.add_argument("-o", "--out_type", default="csv", help="set output file type; default is 'csv' or comma-separated values, a format that can be used by a variety of programs; other valid options are currently 'tsv' or tab-separated values")
    parser.add_argument("-m","--min_count", type=int, default=1, help='minimum count of sequences searched (defaults to 1, keeping all sequences present in the "k-seq start" round, but can be set higher')
    parser.add_argument("-p","--track_progress", action='store_true', help='if this flag is enabled, terminal will output updates on how much progress the code has made')

    args=parser.parse_args()
    
    searchSet = {} #if this is left blank, it will read all sequences at each round; otherwise, will be limited as follows:
        
    with open(args.kseq_rounds) as f:
        roundList = []
        for lineRead in f:
            roundList.append(lineRead.split('\n')[0])
                           
    with open(args.normalization_list) as f:
        normList = []
        for lineRead in f:
            normList.append(float(lineRead))
 
    center = ''
    
    if args.search_set[0] == 'center':
        center = args.search_set[1]
        #if we're only interested in kseq values over one sequence, we limit our search:
        searchSet = readSinglePeak(args.start_round, args.search_set[1], int(args.search_set[2]), args.min_count)

    elif args.search_set[0] == 'list':
        #if we're only interested in kseq values for a specific list of sequences, we limit our search:
        searchSet = readSeqList(args.search_set[1])
        
    (counts, uniqs, tots) = readAll(args.start_round, roundList, normList, args.min_count, 1, searchSet, args.in_type, args.track_progress)
    #read all sequence abundances in all rounds
    
    with open(args.substrate_concs) as f:
        substList = []
        for lineRead in f:
            substList.append(float(lineRead))

    with open(args.rounds_to_average) as f:
        avgList = []
        for lineRead in f:
            avgSublist = []
            
            elems = lineRead.split(',')
            for elem in elems:
                avgSublist.append(int(elem))
            avgList.append(avgSublist)

    with open(args.rounds_to_error) as f:
        errList = []
        for lineRead in f:
            errSublist = []
            
            elems = lineRead.split(',')
            for elem in elems:
                errSublist.append(int(elem))
            errList.append(errSublist)
    
    separator = ','
    if args.out_type == 'tsv':
        separator = '\t'
    
    calculateKinetics(args.output, counts, uniqs, tots, args.start_round, roundList, substList, separator, avgList, errList, center, args.verbose, args.track_progress)

#Reads a list of sequences, each taking up one row in a text file:
def readSeqList(loc):
    seqSet = {}
    
    with open(loc) as f:
        
        for lineRead in f:
            seqSet[lineRead.split('\n')[0]]=0
            
    return seqSet

#Reads a single peak, from a "sequence counts file." Will look at every sequence within a certain 
def readSinglePeak(loc, center, max_dist, minCount=1, trackProgress=False):
    #fileType: 'counts' refers to a list of sequences followed by count numbers, with three lines of header information
    #minCount is the minimum count threshold accepted
    #max_dist is the maximum allowed distance from a center sequence

    seqSet = {}
    
    with open(loc) as f:
        line0 = next(f)
        next(f)
        next(f)            
        
        uniqueSplit = [elem for elem in line0.strip().split()]
        uniques = int(uniqueSplit[-1]) #total number of unique sequences present in this round
        #for "counts" files, we skip three lines of header data
        
        i=0
        
        for lineRead in f:
            line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
            seq = line[0]

            if int(line[1]) >= minCount: #is the count high enough?
                if Levenshtein.distance(center, seq) <= max_dist: #is the distance low enough?
                        seqSet[seq]=0
                            
            if trackProgress:
                i+=1
                if i%100000 == 0:
                    print "Read " + str(i) + " of " + str(uniques) +" sequences from " + loc
            
    return seqSet

#Reads all sequences in a specific k-seq round, from a "sequence counts file."
#Will look at every sequence, unless given a whiteList dictionary object of sequences to search over 
def readSeqs(loc, fileType='counts', minCount=1, norm2=1, whiteList={}, trackProgress=False):
    #fileType: 'counts' refers to a list of sequences followed by count numbers, with three lines of header information
    #minCount: minimum count in the k-seq start round (anything lower is discarded)
    #norm2: normalization constant for this round (in addition to total number of sequences)
    #whiteList: if we're only interested in a certain subset of sequences, providing them here as a dictionary will limit the search

    allSeqCounts = {}
    with open(loc) as f:
        if fileType == 'counts':
            # get unique and total number of reads from the counts file metadata
            line0 = next(f)
            uniqueSplit = [elem for elem in line0.strip().split()]
            uniques = float(uniqueSplit[-1]) #number of unique sequences present in this round
            
            line1 = next(f)
            totalSplit = [elem for elem in line1.strip().split()]
            totals = float(totalSplit[-1]) #total number of unique sequences present in this round
            
            next(f) #skip a blank line in the counts file
            
            i=0
            
            for lineRead in f:
                
                if trackProgress:
                    i+=1
                    if i%100000 == 0:
                        print "Read " + str(i) + " of " + str(int(uniques)) +" sequences from " + loc

                line = [elem for elem in lineRead.strip().split()]
                if int(line[1]) >= minCount: #is the count high enough?
                    if whiteList: #is this sequence wanted?
                        if line[0] in whiteList:
                            allSeqCounts[line[0]] = float(line[1])/totals/norm2
                    else:
                        allSeqCounts[line[0]] = float(line[1])/totals/norm2 #we're adding counts normalized by both total # of sequences and by an arbitrary second normalization constant 

            return (allSeqCounts, uniques, totals)
            #a list of all sequences we want to search over, and their appearance

#Align all sequences found in the new round and their abundance to those found in previous rounds
def alignCounts(masterCounts, roundCounts, rnd, maxRnds, initialize=False):
#masterCounts: A list of all sequences we care about (any not in this key set will be discarded), along with a list of abundance data for all k-seq rounds we've searched so far
#roundCounts: A dictionary of all sequences in the new round we're adding to masterCounts
#rnd: the round number
#maxRnds: If we're initializing, we need to create a list of this many rounds
    
    if initialize:
        for seq in roundCounts.keys():
            counts = [0]*maxRnds
            counts[0] += roundCounts[seq]
            masterCounts.append([seq, counts]) #masterCounts is a list of tuples of the format (sequence, list of abundances)
    else:
        for seq in masterCounts.keys():
            if seq[0] in roundCounts:
                seq[1][rnd] += roundCounts[seq[0]]

    return masterCounts

#The function that collects data over all kseq rounds
def readAll(firstLoc, otherLoc, normList, firstMin=1, otherMin=1, whiteList={}, fileType="counts", trackProgress=False):
    #firstLoc: The "start" round's file location (counts file)
    #otherLoc: A list of all other rounds' file locations (counts files)
    #firstMin, otherMin: minimum counts allowed for each of these rounds
    #normList: a list of normalization constants for all k-seq rounds (based on expected fraction of sequences that survive in each k-seq round from the initial pool), with the first entry as normalization for the start round
        
    masterCounts = []
    maxRnds = len(otherLoc) + 1 #the number of rounds necessary for storing counts data on
    
    (round0counts, round0uniques, round0totals) = readSeqs(firstLoc, fileType, minCount=firstMin, norm2=normList[0], whiteList=whiteList, trackProgress=trackProgress)
    
    masterUniques = [round0uniques]
    masterTotals = [round0totals]
    masterCounts = alignCounts(masterCounts, round0counts, 0, maxRnds, initialize=True)
    
    masterSet = set([seq[0] for seq in masterCounts])
    print("%i sequences found in input pool" %len(round0counts))

    thisRnd = 0
    for loc in otherLoc:

        thisRnd += 1
        (seqs, uniqs, tots) = readSeqs(loc, fileType, otherMin, normList[thisRnd], whiteList=masterSet, trackProgress=trackProgress)
        
        masterUniques.append(uniqs)
        masterTotals.append(tots)
        masterCounts = alignCounts(masterCounts, seqs, thisRnd, maxRnds)
    
    return (masterCounts, masterUniques, masterTotals)

#This function takes all collected data, and calculates k*t and A kinetic values for all sequences, following [surviving fraction] = A(1-Exp(-k*[S]*t))
def calculateKinetics(outLoc, counts, uniqs, tots, firstLoc, otherLoc, substrateConc, separator=',', rndsToAvg=[], rndsToError=[], centerSeq='', verbose=True, trackProgress=False):
    #outLoc: location/name of output file
    #counts, uniques, tots: output from a previous readAll
    #firstLoc/otherLoc: names of rounds read (to use as column headers)
    #separator: character used to separate columns in output
    #substrateConc: list of substrate concentrations for each set of rounds (sets of rounds for each concentration defined in fitAvg)
    #rndsToAvg: list of sets of rounds to average together for fits (e.g. [[1, 2, 3],[4,5,6]] will average the abundances of rounds 1,2, and 3, and then also average 4,5, and 6). If only a single replicate was carried out, then [[1],[2],[3]] etc. If you just want abundance without kseq calculations, leave blank
    #rndsToError: Sets of rounds used for each replicate (e.g.[[1,3,5],[2,4,6]] will evaluate k-seq for just each of the two separate sets, then use this to find st.dev for k/A). If no replicates, leave blank (and st.dev will be based on goodness of fit instead)
    #centerSeq: if provided, each sequence's data will include a column listing its distance to the center sequence of an investigated peak 
    #verbose: if set to true, will provide all round abundance data as well as kseq data for all sequences; if false, will only provide name and kseq data 
    
    #Setting up some header data:
    line0 = 'Seq. Name'
    if verbose:
        line0 += (separator + 'Sequence amount ' + firstLoc)
        for loc in otherLoc:
            line0 += (separator + 'Surv. fraction ' + loc)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    if rndsToAvg:
        line0 += (separator + 'A by avg' + separator + 'k*t by avg')
        line0 += (separator + 'A st.dev' + separator + 'k*t st.dev')
    if centerSeq:
        line0 += (separator + 'distance from ctr')
    line0 += '\n'
    
    #Setting up some round info:
    line1 = ('Unique Seq.s')
    if verbose:
        line1 += (separator + str(uniqs[0]))
        for uniq in uniqs:
            line1 += (separator + str(uniq))
    line1 += '\n'
        
    line2 = ('Total Seq.s')
    if verbose:
        line2 += (separator + str(tots[0]))
        for tot in tots:
            line2 += (separator + str(tot))
    line2 += '\n'

    if rndsToAvg:
        xdata = np.array(substrateConc)
        
        if len(substrateConc) != len(rndsToAvg):
            print "concentrations vs. rounds to average length mismatch error"
    
    with open(outLoc,'w') as outFile:
        outFile.write(line0)
        outFile.write(line1)
        outFile.write(line2)
        
        def func(x, A, k):
            return A * (1 - np.exp(-k * x) )
    
        total_seqs = len(counts)
        seqs_counted = 0
                
        for seq in counts:
            
            if trackProgress:
                seqs_counted += 1
                if seqs_counted % 1000 == 0:
                    print ("calculated " + str(seqs_counted) + " out of " + str(int(total_seqs)) + " sequences")
            
            lineOut = seq[0] #seq's sequence value
            
            if verbose:
                lineOut += (separator + str(seq[1][0]))
                for value in seq[1][1:]: #print surviving fraction of sequences for each k-seq round, normalized by start round abundance
                    if seq[1][0]>0:
                        lineOut += (separator + str(value/(float(seq[1][0]))))
                    else:
                        lineOut += (separator + str(value))
                
            if rndsToAvg:
                avgs = []
                for timePoint in rndsToAvg:
                    timePtCounts = []
                    for rnd in timePoint:
                        if seq[1][0]>0:
                            timePtCounts.append(seq[1][rnd]/(float(seq[1][0]))) #average all rounds' survival concentrations (treating zeros as 0-valued)
                        else:
                            timePtCounts.append(1)
                    avgs.append(np.average(timePtCounts))
                ydata = np.array(avgs)
                
                try:
                    popt, pcov = curve_fit(func, xdata, ydata, method='trf', bounds=(0, [1., np.inf])) #fit k and A from xdata (concentration array) and ydata (reacted fraction array)
                except RuntimeError:
                    popt = [0,0]
                
                lineOut += (separator + str(popt[0]) + separator + str(popt[1]))
                                
                Aset = []
                kset = []

                if popt[1] > 0.0001:
                    #arbitrary threshold? but if popt[1] is too low it's effectively 0
                
                    if rndsToError:      
                        for rndSet in rndsToError:
                            setVals = []
                            for rnd in rndSet:
                                if seq[1][0]>0:
                                    setVals.append(seq[1][rnd]/float(seq[1][0]))
                                else:
                                    setVals.append(1)
                            popt, pcov = curve_fit(func, xdata, np.array(setVals), method='trf', bounds=(0, [1., np.inf]))  #curve fit each set separately
                            Aset.append(popt[0])
                            kset.append(popt[1])
                        
                        lineOut += (separator + str(np.std(Aset)) + separator + str(np.std(kset))) #take the st.dev of all k and A values calculated 
                      
                    else: #if only a single replicate is used, error in k/A is based on fitting parameters instead
                        perr = np.sqrt(np.diag(pcov))
                        lineOut += (separator + str(perr[0]) + separator + str(perr[1]))
                        
                else:
                    lineOut += (separator + '0' + separator + '0')
                    
            if centerSeq:
                lineOut += (separator + str(Levenshtein.distance(centerSeq,seq[0]))) #print center sequences if necessary
                
            lineOut += '\n'
            
            #print lineOut
            outFile.write(lineOut)
    
if __name__ == "__main__":
    main()
