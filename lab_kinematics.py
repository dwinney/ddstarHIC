# Convert .dat files from simulation to a ROOT file and calculate energies in lab frame
#
# Author:       Daniel Winney (2020)
# Affiliation:  Joint Physics Analysis Center (JPAC)
# Email:        dwinney@iu.edu
# ---------------------------------------------------------------------------

import sys, math
from histograms import make_histograms
from data_file import create_data_file
try:
    import ROOT
except:
    exit()
try:
    import numpy as np
except:
    exit()

from data_file import data_file

def main(argv):
    if len(argv) == 0:
        print('Usage: lab_kinematics.py <datfile> <CM energy in GeV>')
        sys.exit(0)
    elif len(argv) == 1:
        infilename = argv[0]
        s = 2.76E3 
    else:
        infilename = argv[0]
        s = float(argv[1]) 
    
    outfilename = infilename.replace('.dat', '.root')
    file = create_data_file(s, infilename, outfilename)
    # make_histograms(file)
    file.Close()

main(sys.argv[1:])