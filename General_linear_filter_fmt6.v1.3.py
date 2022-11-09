##########################################################
### Import Necessary Modules #############################

import argparse                       #provides options at the command line
import sys                            #take command line arguments and uses it in the script
import gzip                           #allows gzipped files to be read
import re                             #allows regular expressions to be used
import Linear_Alignments_v4          #allows linear alignments to be found and retained
#import Linear_Alignments_v4pid          #allows linear alignments to be found and retained

##########################################################
### Command-line Arguments ###############################

parser = argparse.ArgumentParser(description="A script to remove non-linear alignments.")
parser.add_argument("-aln", help = "The location of the alignment file (blast fmt = 6)", default=sys.stdin, required=True)
parser.add_argument("-gap", help = "The maximum distance a linear alignment can be apart (in basepair of query), default:1000", default=1000)
parser.add_argument("-min", help = "The minimum basepair that a linear alignment can be to be retained (includes overlapping bp), default:0", default=0)
parser.add_argument("-pid", help = "The minimum average percent identity to keep a linear alignment, default:75", default=75)
parser.add_argument("-mpid", help = "The minimum percent identity to keep any partial alignment, default:75", default=75)
parser.add_argument("-alnfmt", help = "The alignment format, default:alnfmt6, option:lastz", default="alnfmt6")
parser.add_argument("-self", help = "Allow alignments if query and subject have the same id, default = no", default="no")
parser.add_argument("-print", help = "Print individual alignments (to std. out) or just print summary regions (to std. err), default = yes", default="yes")
#parser.add_argument("-outFmt", help = "For the region summaries, output the region or the query percent IDs, default = region, option pid", default="region")
args = parser.parse_args()

#########################################################
###Variables ############################################

class Variables():
   scaffoldSizes = {}

#########################################################
### Body of script ######################################

class OpenFile():
    def __init__ (self, f, typ, occ):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r') 
        if typ == "aln":
            sys.stderr.write("\nOpened alignment file: {}\n".format(occ))
            OpenAln(self.filename,occ)

class OpenAln():
    def __init__ (self,f,o):
        """Reads an alignment file produced by BLAST, format 6, and outputs all unique query/subject pairs to another subroutine"""
        self.aln = "NA"
        self.alignments = ""
        self.alignments_formatted = ""
        self.maxLength = 0
        self.lineCount = 0
        self.totalLine = 0
        self.linearAlignmentCount = 0
        self.keptLines = 0
        for self.line in f:
            ### Allows gzipped files to be read ###
            try:
                line = line.decode('utf-8')
            except:
                pass          
            if not re.search("^#", self.line):
                self.line = self.line.rstrip('\n')
                if args.alnfmt == "alnfmt6":
                    self.query, self.subject, self.pi, self.length, self.mismatch, self.gap, self.qStart, self.qEnd, self.sStart, self.sEnd = self.line.split("\t")[0:10]

                #########################
                #### Lastz formating ####
                elif args.alnfmt == "lastz":
                    self.query, self.subject, self.identity, self.pid, self.nMatch, self.coverage, self.covPer, self.ngap, self.qStart, self.qEnd, self.sStart, self.sEnd, self.qOrient, self.sOrient = self.line.split("\t")
                    self.pi = self.pid[:-1]
                    self.length = self.coverage.split("/")[0]
                    self.qLength = int(self.coverage.split("/")[1])
                    if self.sOrient == "-":
                        self.tmp = self.sStart
                        self.sStart = self.sEnd
                        self.sEnd = self.tmp
                    self.line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.query, self.subject, self.pi, 
                    self.length, self.coverage, self.covPer, self.qStart, self.qEnd, self.sStart, self.sEnd, self.qOrient, self.sOrient)
                #########################
                #########################
                
                ### Progress on screen ###
                ### and line counted #####
                self.lineCount += 1
                self.totalLine += 1
                if self.lineCount == 1000000:
                    self.lineCount = 0
                    #sys.stderr.write("Processed {} lines\n".format(self.totalLine))
                
                ##########################    
                ### Starting to filter ###
                ##########################               
                ### go to next line if the query and subject are the same ###
                if args.self == "no" and self.query == self.subject:
                    continue                 
                ### only continue with alignments with high enough %id ###
                if float(self.pi) >= float(args.mpid):
                    ### Add alignment to list of alignments ###
                    if self.aln == "NA":
                        self.alignments = "{}".format(self.line)
                        self.alignments_formatted = "{}\t{}\t{}\t{}\t{}\t{}".format(self.qStart, self.qEnd, self.sStart, self.sEnd, self.length, self.pi)
                    ### If it is the end of this query subject, then send it to the linear subroutine ###
                    elif self.aln != "{}\t{}".format(self.query, self.subject):
                        ### Only send if it has the potential to be longer than the minimum ###
                        if int(self.maxLength) >= int(args.min):
                            self.Linear(self.alignments, self.alignments_formatted, self.aln)
                        self.maxLength = 0
                        self.alignments = "{}".format(self.line)
                        self.alignments_formatted = "{}\t{}\t{}\t{}\t{}\t{}".format(self.qStart, self.qEnd, self.sStart, self.sEnd, self.length, self.pi)
                    else:
                        self.alignments += "\n{}".format(self.line)
                        self.alignments_formatted += ",{}\t{}\t{}\t{}\t{}\t{}".format(self.qStart, self.qEnd, self.sStart, self.sEnd, self.length, self.pi)
                    self.maxLength += int(self.length)
                    self.aln = "{}\t{}".format(self.query, self.subject)
        ### After the last line of the file is read, this would be the last alignment set to be analyzed ###            
        if int(self.maxLength) >= int(args.min):
            self.Linear(self.alignments, self.alignments_formatted, self.aln)
        sys.stderr.write("\tFound {} alignments\n".format(self.totalLine))

    def Linear(self, a, af, queRef):
        """Identifies linear alignments passing the criteria in the command line options and outputs them"""          
        self.alignments_formatted = af
        #if args.outFmt == "region":
        self.retained = Linear_Alignments_v4.linearAlignments(self.alignments_formatted, args.gap, args.min, args.pid, queRef)
        #else:
        #    self.retained = Linear_Alignments_v4pid.linearAlignments(self.alignments_formatted, args.gap, args.min, args.pid, queRef)
        if self.retained != "NA":
            if args.print == "yes":
                for self.keep in self.retained:
                    print (a.split("\n")[int(self.keep)])

### Order of script ####
if __name__ == '__main__':
    Variables()
    open_aln = OpenFile(args.aln, "aln", args.aln)
