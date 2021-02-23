# Convert .dat files from simulation to a ROOT file and calculate energies in lab frame
#
# Author:       Daniel Winney (2020)
# Affiliation:  Joint Physics Analysis Center (JPAC)
# Email:        dwinney@iu.edu
# ---------------------------------------------------------------------------

import sys, os
from data_file import data_file
try:
    import ROOT
except:
    exit()
try:
    import numpy as np
except:
    exit()

def main(argv):
    if len(argv) == 0:
        s = 2.76E3 
    elif len(argv) == 1:
        s = float(argv[0]) 
    else:
        sys.exit(0)
    
    home_dir = os.getcwd()

    # Grab all .dat files in /data/ folder
    infiles = []
    for file in os.listdir(home_dir + "/data"):
        if file.endswith(".dat"):
            infiles.append(home_dir + "/data/" + file)
    
    # Import all different files into a root file 
    outfilename = home_dir + "/ddstarpairs.root"
    datafile = data_file(s, infiles, outfilename)
    datafile.close()

main(sys.argv[1:])