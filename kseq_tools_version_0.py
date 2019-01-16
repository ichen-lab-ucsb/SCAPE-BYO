

#Note: the current script is a placeholder/work in progress, meant to be tinkered with
#during use/analysis. It cannot currently be run from the command line, and 
#contains a number of diagnostic outputs to the console window, as well as
#additional functionality not necessary for main use.
#(An example of use is included, commented out, as an example of use)

#A new version, cleaned up and with better commenting, is being worked on.
#This version is included for use in the meanwhile, as necessary.


import numpy as np
from scipy.optimize import curve_fit
import Levenshtein
import time
from sys import getsizeof


def readRndOverlaps(loc1, loc2, outLoc):

    testVal1 = 0
    testVal2 = 0
    testVal3 = 0
    
    seqDict = {}
    seqList = []
    with open(loc1) as f1:
        next(f1)
        next(f1)
        next(f1)
        
        for lineRead in f1:
            testVal1 += 1
            if testVal1 % 10000 == 0:
                print testVal1
            line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
            seqDict[line[0]] = int(line[1])
            
    with open(loc2) as f2:
        next(f2)
        next(f2)
        next(f2)
        
        for lineRead in f2:
            testVal2 += 1
            if testVal2 % 100000 == 0:
                print testVal2

                
            line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
            if line[0] in seqDict:
                seqList.append(line[0])
                testVal3 += 1
                if testVal3 % 100 == 0:
                    print testVal3
                
    with open(outLoc,'w') as outFile:
        for seq in seqList:
            outFile.write(seq + '\n')


                

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

def readSeqList(loc):
    seqList = []
    
    with open(loc) as f:
        
    #    next(f)
    #    next(f)
    #    next(f)
        
        i = 0
        
        for lineRead in f:
            seqList.append(lineRead.split('\n')[0])
            
            i += 1
            if i%1000 == 0:
                print i
   
    print seqList
        
    return seqList

def readSinglePeak(loc, center, min_dist, minCount=1, read_with_count=False):
    seqList = []
    
    with open(loc) as f:
        next(f)
        next(f)
        next(f)
        
        testval = 0
        
        for lineRead in f:
            line = [elem for elem in (lineRead.split(' ')) if (elem != '' and elem !='\n' )]
            seq = line[0]
            if len(seq) > 152:
                seq = seq[0][4:-4]

            testval += 1
            if int(line[1]) >= minCount:
                if Levenshtein.distance(center, seq) <= min_dist:
                    
                    if read_with_count:
                        seqList.append([seq,str(int(line[1]))])
                    else:
                        seqList.append(seq)
            
            if testval % 10000 == 0:
                print testval
                
    
                
                

    return seqList


def readSeqs(loc, fileType='counts', cutOff=3, normalize=True, norm2=1, whiteList={}, similarDistance=0, fake_test_variable=0):
    #fileType: 'counts' refers to a list of sequences followed by count numbers, with three lines of header information
        
    allSeqCounts = {}

    validSet = {}
    print '00000'

    if fake_test_variable == 1:
        with open('testfile-1.txt', 'w') as outfile:
            outfile.write(str('yesss'))

    
    #print whiteList
    
    if fake_test_variable == 1:
        with open('testfile-2.txt', 'w') as outfile:
            outfile.write(str('yesss'))

            print 'oh my'
    
    for seq in whiteList:
        validSet[seq] = 'dummy'
    

    
    print loc
    
    with open(loc) as f:
        if fileType == 'counts':
            
            
            line0 = next(f)
            uniqueSplit = [elem for elem in line0.strip().split()]
            print uniqueSplit           
            uniques = float(uniqueSplit[-1])
            
            line1 = next(f)
            totalSplit = [elem for elem in line1.strip().split()]
            print totalSplit
            totals = float(totalSplit[-1])
            
            next(f)
            #skip3 lines total
            
            i = 0
            
            
            
            for lineRead in f:
                
                timepoint = time.clock()
                
                line = [elem for elem in lineRead.strip().split()]


                i += 1
                
                if i%100000 == 0:
                    print i
                    print loc
                    print "seqC"
                    print getsizeof(validSet)
                    print len(validSet)
                    print getsizeof(allSeqCounts)                    
                    print len(allSeqCounts)

                if i == 1:
                    print line
                if int(line[1]) >= cutOff:
                    
                    if len(line[0]) > 152:
                        line[0] = line[0][4:-4]
                    
                    if similarDistance>0:
                        foundDist = similarDistance + 1
                        bestSeq = ''
                                           

                        for listSeq in whiteList:
                            if Levenshtein.distance(line[0][0:10],listSeq[0:10]) < similarDistance:
                                if Levenshtein.distance(line[0][0:30],listSeq[0:30]) < similarDistance:
                                    distA = Levenshtein.distance(line[0],listSeq)
                                    if distA < foundDist:
                                        foundDist = distA
                                        bestSeq = listSeq
                                        

                        if foundDist <= similarDistance:
                            if bestSeq in allSeqCounts:
                                allSeqCounts[bestSeq] += float(line[1])/totals/norm2
                            else:
                                allSeqCounts[bestSeq] = float(line[1])/totals/norm2
                                        
                            
                    
                    else:
                        
                        if whiteList:
                            
                            if line[0] in validSet:
                                if normalize:
                                    if   line[0] in allSeqCounts:
                                        allSeqCounts[line[0]] += float(line[1])/totals/norm2                                        
                                    else:
                                        allSeqCounts[line[0]] = float(line[1])/totals/norm2
                                        
                                    if line[0]=='AGAAACACGTTGGAATGATGGGAATGTCCCATTCCGCGGTCTAAATATGAGTACAGAATAAGAGAACGTTTATATAATAGTGCTCTTGAATAATCCGTTACGGCCGGGATCGACCACGAGATGCCATCGGTGGTTCTTCAGTTTGGCCCT':
                                        print 'danger will robinson!'
                                        print lineRead
                                        print line[1]
                                        print allSeqCounts[line[0]]
                                        print allSeqCounts[line[0]]*totals*norm2
                                    
                                else:
                                    allSeqCounts[line[0]] = int(line[1])
        
                        else:
                            if normalize:
                                allSeqCounts[line[0]] = float(line[1])/totals/norm2
                            else:
                                allSeqCounts[line[0]] = int(line[1])
                
                                                                   
                    
    return (allSeqCounts, uniques, totals)
    #a list of all sequences we want to search over, and their appearance 

def alignCounts(masterCounts, roundCounts, rnd, maxRnds, initialize=False):

    if initialize:
        for seq in roundCounts:
            counts = [0]*maxRnds
            counts[0] += roundCounts[seq]
            masterCounts.append([seq, counts])
    else:
        for seq in masterCounts:
            if seq[0] in roundCounts:
                seq[1][rnd] += roundCounts[seq[0]]

    return masterCounts
        
def addOffCounts(masterCounts, loc):

    (counts, un, to) = readSeqs(loc, cutOff=1, normalize=False)

    for seq in masterCounts:
        if seq[0] in counts:
            seq.append(counts[seq[0]])
    return masterCounts


def readAll(firstLoc, otherLoc, firstMin, otherMin, peakCtrs=False, searchList=[], normList=[], similar=0):
    
    
    if peakCtrs:
        searchList.append(readPeakCenters(firstLoc))
    
#    print searchList
    
    print 'qqq0'
    
    masterCounts = []
    maxRnds = len(otherLoc) + 1
    if normList:
        (round0Counts, round0uniques, round0totals) = readSeqs(firstLoc, cutOff=firstMin, normalize=True, norm2 = normList[0], whiteList=searchList, similarDistance=similar)
    else:    
        (round0Counts, round0uniques, round0totals) = readSeqs(firstLoc, cutOff=firstMin, normalize=True, whiteList=searchList, similarDistance=similar)
    masterUniques = [round0uniques]
    masterTotals = [round0totals]
    
    print 'qqq1'
    
    if searchList:
        for seq in searchList:
            masterCounts.append([seq, [0]*maxRnds])
            
    print 'qqq2'
    
    masterCounts = alignCounts(masterCounts, round0Counts, 0, maxRnds, initialize=(searchList==[]))
    masterSet = {}
    for seq in masterCounts:
        masterSet[seq[0]] = 0
        



    print len(round0Counts)
    print len(masterCounts)
    
    print 'neato'
    
    thisRnd = 0
    for loc in otherLoc:
        
        print 'qqq3'
        
        thisRnd += 1
        print thisRnd
        if normList:
            if searchList:
                (seqs, uniqs, tots) = readSeqs(loc, cutOff=otherMin, norm2=normList[thisRnd], whiteList=searchList, similarDistance=similar, fake_test_variable=1)
            else:
                (seqs, uniqs, tots) = readSeqs(loc, cutOff=otherMin, norm2=normList[thisRnd], whiteList=masterSet, similarDistance=similar, fake_test_variable=1)

        else:
            if searchList:
                (seqs, uniqs, tots) = readSeqs(loc, cutOff=otherMin, whiteList=searchList, similarDistance=similar, fake_test_variable=1)
            else:
                (seqs, uniqs, tots) = readSeqs(loc, cutOff=otherMin, similarDistance=similar, fake_test_variable=1)
            
        masterUniques.append(uniqs)
        masterTotals.append(tots)
        masterCounts = alignCounts(masterCounts, seqs, thisRnd, maxRnds)
    
    return (masterCounts, masterUniques, masterTotals)
    
def printCounts(outLoc, counts, uniqs, tots, firstLoc, otherLoc, separator='\t', fitAvg=[], centerSeq='', normByInit=False, compact=False, replaceCounts=False):
    
    
    if compact == 1:
        line0 = firstLoc
    if compact == 2:
        line0 = ('X' + separator + 'Abun ' + firstLoc)
    else:
        line0 = ('X' + separator + firstLoc)
        for loc in otherLoc:
            print firstLoc
            print loc
            line0 += (separator + loc)        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    if fitAvg:
        line0 += (separator + 'L by avg' + separator + 'k by avg' + separator + 'L stdev' + separator + 'k stdev')
    if centerSeq:
        line0 += (separator + 'distance from ctr')
    line0 += '\n'
    
    
    line1 = ('Unique')
    
    if compact == 2:
        line1 += (separator + str(uniqs[0]))
    else:    
        for uniq in uniqs:
            line1 += (separator + str(uniq))
    line1 += '\n'
        
        
    line2 = ('Total')
    if compact == 2:
        line2 += (separator + str(tots[0]))
    else:
        for tot in tots:
            line2 += (separator + str(tot))
    line2 += '\n'

    if fitAvg:
        xdata = np.array(fitAvg[0])
        rndsToAvg = fitAvg[1]
        rndsToError = fitAvg[2]
        
        if len(fitAvg[0]) != len(fitAvg[1]):
            print "length mismatch error"
    
    with open(outLoc,'w') as outFile:
        outFile.write(line0)
        outFile.write(line1)
        outFile.write(line2)
        
        def func(x, L, k):
            return L * (1 - np.exp(-k * x) )
    
        counter_o_counts = 0
        
        print len(counts)
        
        for seq in counts:
            
            
            counter_o_counts += 1
            if counter_o_counts % 1000 == 0:
                print counter_o_counts
                print 'yeaaaah boiii!'
            
            if compact == 1:
                if replaceCounts:
                    lineOut = str(seq[2])
                else:
                    lineOut = str(seq[1][0])
            elif compact == 2:
                lineOut = seq[0]
                lineOut += separator + str(seq[1][0])
            else:    
                #print 'heyyo'
                lineOut = seq[0]
                lineOut += separator + str(seq[1][0])            
                for value in seq[1][1:]:
                #    print value
                    if normByInit:
                        if seq[1][0]>0:
                            lineOut += (separator + str(value/(float(seq[1][0]))))
                        else:
                            lineOut += (separator + str(value))
                    else:
                        lineOut += (separator + str(value))
                #print lineOut
                #print 'compact?' + str(compact)
                
            if fitAvg:
                avgs = []
                for timePoint in rndsToAvg:
                    timePtCounts = []
                    for rnd in timePoint:
                        if seq[1][0]>0:
                            timePtCounts.append(seq[1][rnd]/(float(seq[1][0])))
                        else:
                            timePtCounts.append(1)
                    avgs.append(np.average(timePtCounts))
                ydata = np.array(avgs)
                
                
                try:
                    popt, pcov = curve_fit(func, xdata, ydata, method='trf', bounds=(0, [1., np.inf]))
                except RuntimeError:
                    popt = [0,0]
                
                
                lineOut += (separator + str(popt[0]) + separator + str(popt[1]))
                                
                Lset = []
                kset = []

                #print popt
                if popt[1] > 0.0001:
                    
                    #arbitrary threshold?
                
                    if rndsToError == 'value':
                        
                        perr = np.sqrt(np.diag(pcov))
                        lineOut += (separator + str(perr[0]) + separator + str(perr[1]))
                                                
                    else:
                                            
                        for rndSet in rndsToError:
                            setVals = []
                            for rnd in rndSet:
                                if seq[1][0]>0:
                                    setVals.append(seq[1][rnd]/float(seq[1][0]))
                                else:
                                    setVals.append(1)
                            popt, pcov = curve_fit(func, xdata, np.array(setVals), method='trf', bounds=(0, [1., np.inf]))
                            Lset.append(popt[0])
                            kset.append(popt[1])
                        
                        lineOut += (separator + str(np.std(Lset)) + separator + str(np.std(kset)))
                      
                else:
                    lineOut += (separator + '0' + separator + '0')
                    
#                print seq
#                print Lset
#                print kset
#                print np.std(Lset)
#                print np.std(kset)
                    
                
            if centerSeq:
                lineOut += (separator + str(Levenshtein.distance(centerSeq,seq[0])))
                
            lineOut += '\n'
            
            #print lineOut
            outFile.write(lineOut)
        

#in the test case provided below, a set of data from multiple kseq observations (letters) of different selection conditions (numbered)
#are used, along with normalization data, to build a kseq profile for all sequences from a single peak

#testList = ['counts-1A.txt','counts-1B.txt','counts-1C.txt','counts-1D.txt','counts-1E.txt','counts-1F.txt','counts-2A.txt','counts-2B.txt','counts-2C.txt','counts-2D.txt','counts-2E.txt','counts-2F.txt','counts-3A.txt','counts-3B.txt','counts-3C.txt','counts-3D.txt','counts-3E.txt','counts-3F.txt','counts-4A.txt','counts-4B.txt','counts-4C.txt','counts-4D.txt','counts-4E.txt','counts-4F.txt',]
#testNorm = [0.0005, 0.023823133, 0.023823133, 0.023823133, 0.023823133, 0.023823133, 0.023823133, 0.062784812, 0.062784812, 0.062784812, 0.062784812, 0.062784812, 0.062784812, 0.159915207, 0.159915207, 0.159915207, 0.159915207, 0.159915207, 0.159915207, 0.53032596, 0.53032596, 0.53032596, 0.53032596, 0.53032596, 0.53032596]
#testAvg = ([0.00025, 0.00005, 0.00001, 0.000002], [[1, 2, 3, 4, 5, 6], [7, 8, 10, 11, 12], [13, 14, 15], [19, 20, 21, 23]], [[1, 7, 13, 19], [2, 8, 14, 20], [3, 12, 15, 21]])


#center = 'ATTCACCTAGGTCATCGGGTG'
#testSeqs = readSinglePeak('R5c-counts.txt', center, 3)
#(tempCounts, un, to) = readAll('R5c-counts.txt', testList, 1, 1, searchList=testSeqs, normList=testNorm)
#printCounts ('pk-10-kseqs.csv', tempCounts, un, to, 'R5c-counts.txt', testList, separator=',', fitAvg=testAvg, normByInit=True, centerSeq=center)


