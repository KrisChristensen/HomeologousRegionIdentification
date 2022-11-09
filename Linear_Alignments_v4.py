##########################################################
### Import Necessary Modules

import sys                        #take command line arguments and uses it in the script
import re                         #allows regular expressions to be used
import GeneralOverlap_v1          #Checks for general overlaps

"""
This script is meant to be called from another script.  It retains alignments or anything with the proper format that are linear in their progression.
"""
##################
### Variables ####
sort = {}
##################

##################
### Main body ####
def linearAlignments(input_string, gap, minBP, pid, queRef):
    """Identify start of pathways"""    
    positions = input_string.split(",")
    
    #################################################################
    #################################################################
    ###First put into a dictionary based on the reference start position (this will be used to sort later)
    sort.clear()
    alnCount = 0
    for pos in positions:
        qStart, qEnd, sStart, sEnd, length, perID = pos.split("\t")
        if int(qStart) in sort:
            sort[int(qStart)] += ",{}\t{}\t{}\t{}\t{}\t{}\t{}".format(qStart, qEnd, sStart, sEnd, length, perID, alnCount)
        else:
            sort[int(qStart)] = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(qStart, qEnd, sStart, sEnd, length, perID, alnCount)
        alnCount +=1
        
    ################################################################
    ################################################################
    ###Calls to subroutines
    ###Sort by chromosome position, then start creating linear alignments within a max distance set above
    keepPositions = linear(gap, minBP, pid, queRef)
    if len(keepPositions) > 0:
        return(keepPositions)
    else:
        return("NA")

def linear(g, mbp, pi, queRef):
    """Find all linear pathways""" 
    gap = g
    minBP = mbp
    minPid = pi
    query = queRef.split()[0]
    reference = queRef.split()[1]
    #################################################################
    ### A hash is created with possible linear paths ################
    #################################################################
    currentPath = []
    pathways = {}
    pathwayMetrics = {}
    startNewPath = 1
    pathwayCount = 0
    for scaffoldStart in sorted(sort, key=int):
       #sys.stderr.write("{}\t{}\n".format(scaffoldStart, sort[scaffoldStart]))
       
       ### Check to see if there are multiple alignents at the query start position and choose the longest ###
       if len(sort[scaffoldStart].split(",")) > 1:
           sort[scaffoldStart] = longest(sort[scaffoldStart])
       qStart, qEnd, sStart, sEnd, length, perID, alnCount = sort[scaffoldStart].split("\t")
       
       ###########################################################################################################
       ###Add to previous path or decide to create a new path ####################################################
       ###Note: subject changes orientation with blast
       if startNewPath == 0:
           queryPathStart = int(currentPath[0]) - int(gap)
           queryPathEnd = int(currentPath[1]) + int(gap)
           subjectPathStart = int(currentPath[2]) - int(gap)
           subjectPathEnd = int(currentPath[3]) + int(gap)
           #sys.stderr.write("\t\t\t\t\t{}\t{}\n".format(subjectPathStart, subjectPathEnd))
           if int(currentPath[0]) > int(currentPath[1]) and int(qStart) > int(qEnd):
               queryPathStart = int(currentPath[1]) - int(gap)
               queryPathEnd = int(currentPath[0]) + int(gap)
           if int(currentPath[2]) > int(currentPath[3]) and int(sStart) > int(sEnd):
               subjectPathStart = int(currentPath[3]) - int(gap)
               subjectPathEnd = int(currentPath[2]) + int(gap)
               #sys.stderr.write("\t\t\t\t\t{}\t{}\t{}\t{}\t{}\t{}\n".format(currentPath[2], currentPath[3], sStart, sEnd, subjectPathStart, subjectPathEnd))
           qOverlapType, qOverlapAmount = GeneralOverlap_v1.Overlap(qStart, qEnd, queryPathStart, queryPathEnd)
           sOverlapType, sOverlapAmount = GeneralOverlap_v1.Overlap(sStart, sEnd, subjectPathStart, subjectPathEnd)                      
           #sys.stderr.write("\t{}\n\t{}\n".format(currentPath, sort[scaffoldStart].split("\t")[0:4]))
           #sys.stderr.write("\t\t{}\t{}\t{}\t{}\n".format(queryPathStart, queryPathEnd, subjectPathStart, subjectPathEnd))
           #sys.stderr.write("\t\t\t{}\t{}\t{}\t{}\n".format(qOverlapType, qOverlapAmount, sOverlapType, sOverlapAmount))
           if qOverlapType != "NA" and sOverlapType != "NA":
               pathways[pathwayCount] += ",{}".format(sort[scaffoldStart])
               pathwayMetrics[pathwayCount] = "{}\t{}\t{}".format(int(pathwayMetrics[pathwayCount].split()[0]) + 1, int(pathwayMetrics[pathwayCount].split()[1]) + int(length), float(pathwayMetrics[pathwayCount].split()[2]) + float(perID))
               #sys.stderr.write("{}\n\t\t{}\n".format(pathways[pathwayCount], pathwayMetrics[pathwayCount]))
               #sys.stderr.write("{}\n\t\t{}\n".format(pathwayCount, pathwayMetrics[pathwayCount]))
           else:
               startNewPath = 1
       currentPath = [qStart, qEnd, sStart, sEnd]
       ###########################################################################################################
       ###Check to see if a new path needs to be created #########################################################
       if startNewPath == 1:
           pathwayCount += 1
           pathways[pathwayCount] = "{}".format(sort[scaffoldStart])
           pathwayMetrics[pathwayCount] = "{}\t{}\t{}".format(1,length,perID)
           startNewPath = 0
           continue
    ##### Check metrics before returning ####
    keeps = []
    for pathwayCount in pathwayMetrics:
        totCount = pathwayMetrics[pathwayCount].split()[0]
        totLength = pathwayMetrics[pathwayCount].split()[1]
        avgLength = int(totLength)/float(totCount)
        avgPID = float(pathwayMetrics[pathwayCount].split()[2])/int(totCount)
        if float(avgPID) >= float(minPid) and int(totLength) >= int(minBP):
            #sys.stderr.write("{}\t{}\t{}\t{}\t{}\n".format(pathwayCount, totCount, totLength, avgLength, avgPID))
            #sys.stderr.write("\t{}\n".format(pathways[pathwayCount]))
            startPathway = "NA"
            endPathway = "NA"
            startqPathway = "NA"
            endqPathway = "NA"
            sOrientation = 0
            qOrientation = 0
            ####Section to get metrics about the pathway (start, finish, and orientation)
            for keep in pathways[pathwayCount].split(","):
                qStart, qEnd, sStart, sEnd, length, perID, alnCount = keep.split()
                keeps.append(alnCount)
                tmpStart = sStart
                tmpEnd = sEnd
                tmpqStart = qStart
                tmpqEnd = qEnd
                sOrientation += 1
                qOrientation += 1
                if int(sStart) > int(sEnd):
                    tmpStart = sEnd
                    tmpEnd = sStart
                    sOrientation -= 2
                if int(qStart) > int(qEnd):
                    tmpqStart = qEnd
                    tmpqEnd = qStart
                    qOrientation -= 2
                if startPathway == "NA" or int(tmpStart) < int(startPathway):
                    startPathway = int(tmpStart)
                if endPathway == "NA" or int(tmpEnd) > int(endPathway):
                    endPathway = int(tmpEnd)
                if startqPathway == "NA" or int(tmpqStart) < int(startqPathway):
                    startqPathway = int(tmpqStart)
                if endqPathway == "NA" or int(tmpqEnd) > int(endqPathway):
                    endqPathway = int(tmpqEnd)
            if int(sOrientation) < 0:
                tStart = endPathway
                endPathway = startPathway
                startPathway = tStart
            if int(qOrientation) < 0:
                tStart = endqPathway
                endqPathway = startqPathway
                startqPathway = tStart            
            sys.stderr.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query, startqPathway, endqPathway, reference, startPathway, endPathway, avgPID))
    return(keeps)
            
                
def longest(lon):
    """Return the longest alignment for a specific reference starting position"""        
    bestLen = 0
    bestLenInfo = "NA"
    for aln in lon.split(","):
        qStart, qEnd, sStart, sEnd, length, perID, alnCount = aln.split("\t")
        if int(length) > int(bestLen):
            bestLenInfo = aln
            bestLen = int(length)
    #sys.stderr.write("Original: {}\n\tReturning: {}\n".format(lon, bestLenInfo))
    return(bestLenInfo)
