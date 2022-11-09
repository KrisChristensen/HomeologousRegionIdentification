##########################################################
### Import Necessary Modules #############################

import argparse                       #provides options at the command line
import sys                            #take command line arguments and uses it in the script
import gzip                           #allows gzipped files to be read
import re                             #allows regular expressions to be used

##########################################################
### Command-line Arguments ###############################

parser = argparse.ArgumentParser(description="A script to output the standard output from the General_linear_filter_fmt6.v1.3 script to circos output.")
parser.add_argument("-input", help = "The location of the summary information (in a file) output from the General_linear_filter_fmt6.v1.3 script", default=sys.stdin, required=True)
parser.add_argument("-tbl", help = "The location of a table file that replaces chromosome names and adds colors format: previous chr1<tab>new chr<tab>color", default=sys.stdin, required=True)
parser.add_argument("-tLen", help = "The minimum total length of a region before reporting, default=50000", default=50000)
args = parser.parse_args()

#########################################################
###Variables ############################################

class Variables():
   tbl = {}

#########################################################
### Body of script ######################################

class OpenFile():
    def __init__ (self, f, typ, occ):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r') 
        if typ == "tbl":
            #sys.stderr.write("\nOpened table file: {}\n".format(occ))
            OpenTbl(self.filename,occ)
        elif typ == "input":
            #sys.stderr.write("\nOpened input file: {}\n".format(occ))
            OpenInput(self.filename,occ)

class OpenTbl():
    def __init__ (self,f,o):
        """Reads the table input file with the old chromosome names, the new names, and the colors (tab-delimited)"""
        for self.line in f:
            ### Allows gzipped files to be read ###
            try:
                line = line.decode('utf-8')
            except:
                pass          
            if not re.search("^#", self.line):
                self.line = self.line.rstrip('\n')
                self.oldChr, self.newChr, self.color = self.line.split()
                Variables.tbl[self.oldChr] = "{}\t{}".format(self.newChr, self.color)
        f.close()
                
class OpenInput():
    def __init__ (self,f,o):
        """Reads the input file output from the standard output of the script General_linear_filter_fmt.v1.3 script"""
        self.color = {}
        self.toPrint = {}
        self.toPrintPid = {}
        self.toPrintMetrics = {}
        for self.line in f:
            ### Allows gzipped files to be read ###
            try:
                line = line.decode('utf-8')
            except:
                pass
            if not re.search("^#", self.line):
                self.line = self.line.rstrip('\n')
                try:
                    self.query, self.qStart, self.qEnd, self.reference, self.rStart, self.rEnd, self.pid = self.line.split()
                    if self.query in Variables.tbl and self.reference in Variables.tbl:
                        self.newQuery, self.qColor = Variables.tbl[self.query].split()
                        self.newReference, self.rColor = Variables.tbl[self.reference].split()
                        self.colorToUse = self.qColor
                        if "{}\t{}".format(self.newReference, self.newQuery) in self.color:
                            self.colorToUse = self.color["{}\t{}".format(self.newReference, self.newQuery)]
                        else:
                            self.color["{}\t{}".format(self.newQuery, self.newReference)] = self.qColor
                        #print ("{}\t{}\t{}\t{}\t{}\t{}\tcolor={}".format(self.newQuery, self.qStart, self.qEnd, self.newReference, self.rStart, self.rEnd, self.colorToUse))
                        #sys.stderr.write("{}\t{}\t{}\t{}\n".format(self.newQuery, self.qStart, self.qEnd, self.pid))
                        self.avgLen = int(self.qEnd) - int(self.qStart)
                        self.key = "{}\t{}".format(self.newQuery, self.newReference)                        
                        if self.key in self.toPrintMetrics:
                            self.toPrintMetrics[self.key] += float(self.avgLen)
                            self.toPrint[self.key] += ",{}\t{}\t{}\t{}\t{}\t{}\tcolor={}".format(self.newQuery, self.qStart, self.qEnd, self.newReference, self.rStart, self.rEnd, self.colorToUse)
                            self.toPrintPid[self.key] += ",{}\t{}\t{}\t{}\n".format(self.newQuery, self.qStart, self.qEnd, self.pid)                            
                        else:
                            self.toPrintMetrics[self.key] = float(self.avgLen)
                            self.toPrint[self.key] = "{}\t{}\t{}\t{}\t{}\t{}\tcolor={}".format(self.newQuery, self.qStart, self.qEnd, self.newReference, self.rStart, self.rEnd, self.colorToUse)
                            self.toPrintPid[self.key] = "{}\t{}\t{}\t{}\n".format(self.newQuery, self.qStart, self.qEnd, self.pid)                            
                except:
                    pass
        for self.key1 in self.toPrintMetrics:
            #sys.stderr.write("\t{}\n".format(self.toPrintMetrics[self.key1]))
            if float(self.toPrintMetrics[self.key1]) >= float(args.tLen):
                for self.lineToPrint in self.toPrint[self.key1].split(","):
                    print("{}".format(self.lineToPrint))
                for self.lineToPrint1 in self.toPrintPid[self.key1].split(","):
                    sys.stderr.write("{}".format(self.lineToPrint1))
        f.close()                              

### Order of script ####
if __name__ == '__main__':
    Variables()
    open_tbl = OpenFile(args.tbl, "tbl", args.tbl)    
    open_input = OpenFile(args.input, "input", args.input)
