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

class data_file:
    # constructor 
    # intialize all the data that will be filled into the root tree later
    def __init__(self, _collision_energy, _datfilenames, _outfilename):
        self.s = _collision_energy
        self.outfile = ROOT.TFile( _outfilename, 'RECREATE')

        # Lab frame 
        self.D_pos        =  np.empty((4,1), dtype='float32')
        self.D_mom        =  np.empty((4,1), dtype='float32') 
        self.Dstar_pos    =  np.empty((4,1), dtype='float32') 
        self.Dstar_mom    =  np.empty((4,1), dtype='float32') 
        self.Miss_pos     =  np.empty((4,1), dtype='float32') 
        self.Miss_mom     =  np.empty((4,1), dtype='float32') 

        # CM frame
        self.cm_masses    = np.empty((2,1), dtype='float32')
        self.cm_energies  = np.empty((2,1), dtype='float32')
        self.cm_mom       = np.empty((3,1), dtype='float32')
        self.cm_D_pos     = np.empty((4,1), dtype='float32')
        self.cm_Dstar_pos = np.empty((4,1), dtype='float32')
        
        self.import_data(_datfilenames)
        self.process_data()
    
    # These branches come from the inital raw data file so they're the same everywhere
    def setup_base_branches(self, tree):
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
    def import_data(self, listfiles):
        # ROOT structures
        tree_raw = ROOT.TTree("raw_data", "raw_data")
        self.setup_base_branches(tree_raw)

        for file in listfiles:
            datfile = open(file)
            lines = datfile.readlines()
            for line in lines:

                # convert ascii to doubles
                line_processed = self.process_line(line)
                D, Dstar = line_processed[0], line_processed[1]

                # save values
                self.D_pos[:, 0]     = [0, D[0],     D[1],     D[2]]
                self.D_mom[:, 0]     = [0, D[3],     D[4],     D[5]]
                self.Dstar_pos[:, 0] = [0, Dstar[0], Dstar[1], Dstar[2]]
                self.Dstar_mom[:, 0] = [0, Dstar[3], Dstar[4], Dstar[5]]

                tree_raw.Fill()
            datfile.close()

        tree_raw.Write()
        return

    # Parse each line seperating of the .dat file
    def process_line(self, line):
        line_split = line.split(',')

        D     = [line_split[0], line_split[1], line_split[2], line_split[3], line_split[4],  line_split[5]]
        Dstar = [line_split[6], line_split[7], line_split[8], line_split[9], line_split[10], line_split[11]]

        # Convert to floats
        for i in range(len(D)):
            D[i]     = float(D[i])
            Dstar[i] = float(Dstar[i])

        return [D, Dstar]

    # Take imported data and calculate additional lab quantities:
    # Energies for D and Dstar
    # Assumed t = 0 for lab (will be necessary for boosts)
    # Missing particle 4-position and 4-momentum
    def process_data(self):

        # Lab frame tree
        tree_lab = ROOT.TTree("lab", "lab")
        self.setup_base_branches(tree_lab)
        tree_lab.Branch('D_t',     self.D_pos[0],         "D_t/F")
        tree_lab.Branch('Dstar_t', self.Dstar_pos[0],     "Dstar_t/F")
        tree_lab.Branch('D_E',     self.D_mom[0],         "D_E/F")
        tree_lab.Branch('Dstar_E', self.Dstar_mom[0],     "Dstar_E/F")
        tree_lab.Branch('Miss_E',  self.Miss_mom[0],      "Miss_E/F")
        tree_lab.Branch('Miss_px', self.Miss_mom[1],      "Miss_px/F")
        tree_lab.Branch('Miss_py', self.Miss_mom[2],      "Miss_py/F")
        tree_lab.Branch('Miss_pz', self.Miss_mom[3],      "Miss_pz/F")

        # CM frame tree
        tree_cm = ROOT.TTree("cm", "cm")
        tree_cm.Branch('W_DDstar',  self.cm_masses[0],    "W_DDstar/F")
        tree_cm.Branch('Miss_mass', self.cm_masses[1],    "Miss_mass/F")
        tree_cm.Branch('D_en',      self.cm_energies[0],  "D_en/F")
        tree_cm.Branch('Dstar_en',  self.cm_energies[1],  "Dstar_en/F")
        tree_cm.Branch('CM_px',     self.cm_mom[0],       "CM_px/F")
        tree_cm.Branch('CM_py',     self.cm_mom[1],       "CM_py/F")
        tree_cm.Branch('CM_pz',     self.cm_mom[2],       "CM_pz/F")
        tree_cm.Branch('D_t',       self.cm_D_pos[0],     "D_t/F")
        tree_cm.Branch('D_x',       self.cm_D_pos[1],     "D_x/F")
        tree_cm.Branch('D_y',       self.cm_D_pos[2],     "D_y/F")
        tree_cm.Branch('D_z',       self.cm_D_pos[3],     "D_z/F")
        tree_cm.Branch('Dstar_t',   self.cm_Dstar_pos[0], "Dstar_t")
        tree_cm.Branch('Dstar_x',   self.cm_Dstar_pos[1], "Dstar_x")
        tree_cm.Branch('Dstar_y',   self.cm_Dstar_pos[2], "Dstar_y")
        tree_cm.Branch('Dstar_z',   self.cm_Dstar_pos[3], "Dstar_z")

        # Main loop for each event
        for entry in self.outfile.raw_data:
            self.process_entry(entry)
            tree_lab.Fill()
            tree_cm.Fill()

        tree_lab.Write()
        tree_cm.Write()
        
        return

    def process_entry(self, entry):

        ## 4-positions
        rD_lab, rDstar_lab = ROOT.Math.XYZTVector(), ROOT.Math.XYZTVector()
        rD_lab.SetXYZT(entry.D_x,     entry.D_y,     entry.D_z,     0.)
        rDstar_lab.SetXYZT(entry.Dstar_x, entry.Dstar_y, entry.Dstar_z, 0.)

        ## 4-momenta        
        pD_lab      = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(entry.D_px,     entry.D_py,     entry.D_pz,     mD)
        pDstar_lab  = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(entry.Dstar_px, entry.Dstar_py, entry.Dstar_pz, mDstar)

        ## Missing particle momentum and energy
        pDDstar_lab = (pD_lab + pDstar_lab)
        eMiss_lab   = math.sqrt(self.s) - pD_lab.E() - pDstar_lab.E()
        pMiss_lab   = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(- pDDstar_lab.Px(), - pDDstar_lab.Py(), - pDDstar_lab.Pz(), eMiss_lab)

        # Save lab values
        self.D_pos[:, 0]      = [ rD_lab.T(),     rD_lab.X(),      rD_lab.Y(),      rD_lab.Z()      ]
        self.D_mom[:, 0]      = [ pD_lab.E(),     pD_lab.Px(),     pD_lab.Py(),     pD_lab.Pz()     ]
        self.Dstar_pos[:, 0]  = [ rDstar_lab.T(), rDstar_lab.X(),  rDstar_lab.Y(),  rDstar_lab.Z()  ]
        self.Dstar_mom[:, 0]  = [ pDstar_lab.E(), pDstar_lab.Px(), pDstar_lab.Py(), pDstar_lab.Pz() ]
        self.Miss_mom[:, 0]   = [ pMiss_lab.E(),  pMiss_lab.Px(),  pMiss_lab.Py(),  pMiss_lab.Pz()  ]

        # Now we everything to the CM frame
        boost_vector = pDDstar_lab.BoostToCM()
        rD_cm        = ROOT.Math.VectorUtil.boost(rD_lab,      boost_vector)
        pD_cm        = ROOT.Math.VectorUtil.boost(pD_lab,      boost_vector)
        rDstar_cm    = ROOT.Math.VectorUtil.boost(rDstar_lab,  boost_vector)
        pDstar_cm    = ROOT.Math.VectorUtil.boost(pDstar_lab,  boost_vector)
        pDDstar_cm   = ROOT.Math.VectorUtil.boost(pDDstar_lab, boost_vector)
        pMiss_cm     = ROOT.Math.VectorUtil.boost(pMiss_lab,   boost_vector)

        # Save CM values
        self.cm_energies[:,0]    = [ pD_cm.E(),      pDstar_cm.E() ]
        self.cm_masses[:,0]      = [ pDDstar_cm.E(), pMiss_cm.E()  ]
        self.cm_mom[:,0]         = [ pD_cm.Px(),     pD_cm.Py(),      pD_cm.Pz() ]
        self.cm_D_pos[:,0]       = [ rD_cm.T(),      rD_cm.X(),       rD_cm.Y(),      rD_cm.Z()     ]
        self.cm_Dstar_pos[:,0]   = [ rDstar_cm.T(),  rDstar_cm.X(),   rDstar_cm.Y(),  rDstar_cm.Z() ]

        return