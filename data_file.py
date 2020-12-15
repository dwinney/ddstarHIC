# Container class to sort all the data in the root files
#
# Author:       Daniel Winney (2020)
# Affiliation:  Joint Physics Analysis Center (JPAC)
# Email:        dwinney@iu.edu
# ---------------------------------------------------------------------------

import sys, math
try:
    import ROOT
except:
    exit()
try:
    import numpy as np
except:
    exit()

# PDG Masses in GeV
mD     = 1.86484
mDstar = 2.00685
masses = [mD, mDstar]

def create_data_file(_s, _datfilename, _outfilename):
    d = data_file(_s, _datfilename, _outfilename)
    d.process_data()
    
    return d.outfile

class data_file:

    def __init__(self, _collision_energy, _datfilename, _outfilename):
        self.s = _collision_energy
        self.outfile = ROOT.TFile( _outfilename, 'RECREATE')
        self.D_pos       =  np.empty((4,1), dtype='float32')
        self.D_mom       =  np.empty((4,1), dtype='float32') 
        self.Dstar_pos   =  np.empty((4,1), dtype='float32') 
        self.Dstar_mom   =  np.empty((4,1), dtype='float32') 
        self.Miss_pos    =  np.empty((4,1), dtype='float32') 
        self.Miss_mom    =  np.empty((4,1), dtype='float32') 
        self.import_data(_datfilename)
    
    # These branches come from the inital raw data file so they're the same everywher
    def setup_branches(self, tree):
        tree.Branch('D_x', self.D_pos[1], "D_x/F")
        tree.Branch('D_y', self.D_pos[2], "D_y/F")
        tree.Branch('D_z', self.D_pos[3], "D_z/F")

        tree.Branch('D_px', self.D_mom[1], "D_px/F")
        tree.Branch('D_py', self.D_mom[2], "D_py/F")
        tree.Branch('D_pz', self.D_mom[3], "D_pz/F")

        tree.Branch('Dstar_x', self.Dstar_pos[1], "Dstar_x/F")
        tree.Branch('Dstar_y', self.Dstar_pos[2], "Dstar_y/F")
        tree.Branch('Dstar_z', self.Dstar_pos[3], "Dstar_z/F")

        tree.Branch('Dstar_px', self.Dstar_mom[1], "Dstar_px/F")
        tree.Branch('Dstar_py', self.Dstar_mom[2], "Dstar_py/F")
        tree.Branch('Dstar_pz', self.Dstar_mom[3], "Dstar_pz/F")

    # External option to close file so that data_file can be left open for histogram writing
    def close(self):
        self.outfile.Close()  
        return

    # Parse .dat file and import it to the pre-setup branches
    def import_data(self, datfilename):
        # ROOT structures
        tree_raw = ROOT.TTree("raw_data", "raw_data")
        self.setup_branches(tree_raw)

        datfile = open(datfilename)
        lines = datfile.readlines()
        for line in lines:

            # convert ascii to doubles
            line_processed = self.process_line(line)
            D, Dstar = line_processed[0], line_processed[1]

            # save values
            self.D_pos[:, 0] = [0, D[0], D[1], D[2]]
            self.D_mom[:, 0] = [0, D[3], D[4], D[5]]
            self.Dstar_pos[:, 0] = [0, Dstar[0], Dstar[1], Dstar[2]]
            self.Dstar_mom[:, 0] = [0, Dstar[3], Dstar[4], Dstar[5]]

            tree_raw.Fill()

        tree_raw.Write()
        datfile.close()
        return

    # Parse each line seperating of the .dat file
    def process_line(self, line):
        # some events are D-Dbar others are D-D*bar
        # Here treat them the same but 
        # TODO: filter pairs by particle ID
        line_nolabels = line.replace('421\t-421\t', '')
        line_nolabels = line_nolabels.replace('421\t-423\t', '')

        line_split = line_nolabels.split('}\t{')
        for i, subline in enumerate(line_split):
            subline_new = subline.replace('{', '')
            subline_new = subline_new.replace('}', '')
            subline_new = subline_new.replace('\n', '')
            line_split[i] = subline_new 

        D = line_split[0].split(",")
        for i, entry in enumerate(D):
            D[i] = float(entry)
        Dstar = line_split[1].split(",")
        for i, entry in enumerate(Dstar):
            Dstar[i] = float(entry)

        return [D, Dstar]

    # Take imported data and calculate additional lab quantities:
    # Energies for D and Dstar
    # Assumed t = 0 for lab (will be necessary for boosts)
    # Missing particle 4-position and 4-momentum
    def process_data(self):
        tree_lab = ROOT.TTree("lab", "lab")
        self.setup_branches(tree_lab)

        # New branches for the full 4-vectors
        tree_lab.Branch('D_t',     self.D_pos[0], "D_t/F")
        tree_lab.Branch('Dstar_t', self.Dstar_pos[0], "Dstar_t/F")
        tree_lab.Branch('D_E',     self.D_mom[0], "D_E/F")
        tree_lab.Branch('Dstar_E', self.Dstar_mom[0], "Dstar_E/F")
        tree_lab.Branch('Miss_E',  self.Miss_mom[0], "Miss_E/F")
        tree_lab.Branch('Miss_px', self.Miss_mom[1], "Miss_px/F")
        tree_lab.Branch('Miss_py', self.Miss_mom[2], "Miss_py/F")
        tree_lab.Branch('Miss_pz', self.Miss_mom[3], "Miss_pz/F")

        for entry in self.outfile.raw_data:
            self.process_entry(entry)
            tree_lab.Fill()

        tree_lab.Write()
        return

    def process_entry(self, entry):

        ## Calculate energies for D and Dstar
        # print(entry)
        rD = ROOT.TVector3(entry.D_x, entry.D_y, entry.D_z)
        vpD = ROOT.TVector3(entry.D_px, entry.D_py, entry.D_pz)
        eD  = math.sqrt(vpD.Mag2() + mD**2)

        rDstar = ROOT.TVector3(entry.Dstar_x, entry.Dstar_y, entry.Dstar_z)
        vpDstar = ROOT.TVector3(entry.Dstar_px, entry.Dstar_py, entry.Dstar_pz)
        eDstar = math.sqrt(vpDstar.Mag2() + mDstar**2)

        ## Missing particle momentum and energy
        vpMiss = - (vpD + vpDstar)
        eMiss = math.sqrt(self.s) - eD - eDstar

        # Save values
        self.D_pos[:, 0]        = [ 0, rD.X(), rD.Y(), rD.Z() ]
        self.D_mom[:, 0]        = [eD, vpD.X(), vpD.Y(), vpD.Z()]
        self.Dstar_pos[:, 0]    = [ 0, rDstar.X(), rDstar.Y(), rDstar.Z() ]
        self.Dstar_mom[:, 0]    = [eDstar, vpDstar.X(), vpDstar.Y(), vpDstar.Z()]
        self.Miss_mom[:, 0]     = [eMiss, vpMiss.X(), vpMiss.Y(), vpMiss.Z()]

        return