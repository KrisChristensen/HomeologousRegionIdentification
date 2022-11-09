##########################################################
### Import Necessary Modules

import sys                       #take command line arguments and uses it in the script
import re                       #allows regular expressions to be used

"""
This script is meant to be called from another script.  It identifies if there is an overlap.
"""

def Overlap(s1, e1, s2, e2):      
    status = "NA"
    overlap = 0
    if float(s1) >= float(s2) and float(s1) <= float(e2) and float(e1) >= float(s2) and float(e1) <= float(e2):
        status = "overlap1"
        overlap = float(e1) - float(s1) + 1
    elif float(s2) >= float(s1) and float(s2) <= float(e1) and float(e2) >= float(s1) and float(e2) <= float(e1):
        status = "overlap2"
        overlap = float(e2) - float(s2) + 1
    elif float(s1) > float(s2) and float(s1) < float(e2):
        status = "overlap3"
        overlap = float(e2) - float(s1) + 1
    elif float(s2) > float(s1) and float(s2) < float(e1):
        status = "overlap4"
        overlap = float(e1) - float(s2) + 1
    return (status, overlap)

if __name__ == '__main__':
    start1 = 1
    end1 = 1000
    start2 = 200
    end2 = 2000
    (output, overlap_amount) = Overlap(start1, end1, start2, end2)    
    print ("{}\t{}".format(output, overlap_amount))


####################################################################
####  Overlap1
####
####   s1 ---------- e1
#### s2 ----------------- e2
####
####  Overlap2
####
#### s1 ----------------- e1
####   s2 ---------- e2
####
####  Overlap3
####
####      s1 ------------------- e1
####   s2 ------------- e2
####
####  Overlap4
####
####  s1 ----------- e1
####        s2 ------------- e2
