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
        print("Creating new datafile: " + _outfilename + "\n")
        self.s = _collision_energy
        self.outfile = ROOT.TFile( _outfilename, 'RECREATE')

        # Lab frame 
        self.r1  =  np.empty((4,1), dtype='float32')
        self.r2  =  np.empty((4,1), dtype='float32') 
        self.r3  =  np.empty((4,1), dtype='float32') 
        self.q1  =  np.empty((4,1), dtype='float32') 
        self.q2  =  np.empty((4,1), dtype='float32') 
        self.q3  =  np.empty((4,1), dtype='float32') 

        # CM frame
        self.cm_masses    = np.empty((2,1), dtype='float32')
        self.cm_energy    = np.empty((2,1), dtype='float32')
        self.qcm          = np.empty((3,1), dtype='float32')
        self.rcm1         = np.empty((4,1), dtype='float32')
        self.rcm2         = np.empty((4,1), dtype='float32')
        
        self.import_data(_datfilenames)
        self.process_data()
    
    # These branches come from the inital raw data file so they're the same everywhere
    def setup_base_branches(self, tree):
        tree.Branch('r1_x', self.r1[1], "r1_x/F")
        tree.Branch('r1_y', self.r1[2], "r1_y/F")
        tree.Branch('r1_z', self.r1[3], "r1_z/F")

        tree.Branch('r2_x', self.r2[1], "r2_x/F")
        tree.Branch('r2_y', self.r2[2], "r2_y/F")
        tree.Branch('r2_z', self.r2[3], "r2_z/F")

        tree.Branch('q1_x', self.q1[1], "q1_x/F")
        tree.Branch('q1_y', self.q1[2], "q1_y/F")
        tree.Branch('q1_z', self.q1[3], "q1_z/F")

        tree.Branch('q2_x', self.q2[1], "q2_x/F")
        tree.Branch('q2_y', self.q2[2], "q2_y/F")
        tree.Branch('q2_z', self.q2[3], "q2_z/F")
        
    # External option to close file so that data_file can be left open for histogram writing
    def close(self):
        self.outfile.Close()  
        return

    # Parse .dat file and import it to the pre-setup branches
    def import_data(self, listfiles):
        # ROOT structures
        tree_raw = ROOT.TTree("raw_data", "raw_data")
        self.setup_base_branches(tree_raw)

        print("Importing data from: ")
        nevents = 0
        for file in listfiles:
            print("\t" + file + "...", end =" ")
            datfile = open(file)
            lines = datfile.readlines()
            for line in lines:

                # convert ascii to doubles
                line_processed = self.process_line(line)
                D, Dstar = line_processed[0], line_processed[1]

                # save values
                self.r1[:, 0]     = [0, D[0],     D[1],     D[2]]
                self.q1[:, 0]     = [0, D[3],     D[4],     D[5]]
                self.r2[:, 0] = [0, Dstar[0], Dstar[1], Dstar[2]]
                self.q2[:, 0] = [0, Dstar[3], Dstar[4], Dstar[5]]

                nevents += 1
                tree_raw.Fill()
            print(str(len(lines)) + " events");
            datfile.close()
        tree_raw.Write()

        print("Total: " + str(nevents) + " events imported.")
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
        tree_lab.Branch('r1_t',  self.r1[0],   "r1_t/F")
        tree_lab.Branch('r2_t',  self.r2[0],   "r2_t/F")
        tree_lab.Branch('q1_E',  self.q1[0],   "q1_E/F")
        tree_lab.Branch('q2_E',  self.q2[0],   "q2_E/F")
        tree_lab.Branch('q3_E',  self.q3[0],   "q3_E/F")
        tree_lab.Branch('q3_x',  self.q3[1],   "q3_x/F")
        tree_lab.Branch('q3_y',  self.q3[2],   "q3_y/F")
        tree_lab.Branch('q3_z',  self.q3[3],   "q3_z/F")

        # CM frame tree
        tree_cm = ROOT.TTree("cm", "cm")
        tree_cm.Branch('W12',    self.cm_masses[0],  "W12/F")
        tree_cm.Branch('W3',     self.cm_masses[1],  "W3/F")
        tree_cm.Branch('E1',     self.cm_energy[0],  "E1/F")
        tree_cm.Branch('E2',     self.cm_energy[1],  "E2/F")
        tree_cm.Branch('qcm_x',  self.qcm[0],        "qcm_x/F")
        tree_cm.Branch('qcm_y',  self.qcm[1],        "qcm_y/F")
        tree_cm.Branch('qcm_z',  self.qcm[2],        "qcm_z/F")
        tree_cm.Branch('rcm1_t', self.rcm1[0],       "rcm1_t/F")
        tree_cm.Branch('rcm1_x', self.rcm1[1],       "rcm1_x/F")
        tree_cm.Branch('rcm1_y', self.rcm1[2],       "rcm1_y/F")
        tree_cm.Branch('rcm1_z', self.rcm1[3],       "rcm1_z/F")
        tree_cm.Branch('rcm2_t', self.rcm2[0],       "rcm2_t")
        tree_cm.Branch('rcm2_x', self.rcm2[1],       "rcm2_x")
        tree_cm.Branch('rcm2_y', self.rcm2[2],       "rcm2_y")
        tree_cm.Branch('rcm2_z', self.rcm2[3],       "rcm2_z")

        # Main loop for each event
        print("\nProcessing events:")
        n_processed = 0
        for entry in self.outfile.raw_data:
            n_processed += 1
            self.process_entry(entry)
            tree_lab.Fill()
            tree_cm.Fill()
            if (n_processed % 5000 == 0):
                    print("\t" + str(n_processed) + " events processed...")

        tree_lab.Write()
        tree_cm.Write()
        print("Done!")
        return

    def process_entry(self, entry):

        ## 4-positions
        rD_lab, rDstar_lab = ROOT.Math.XYZTVector(), ROOT.Math.XYZTVector()
        rD_lab.SetXYZT(    entry.r1_x, entry.r1_y, entry.r1_z, 0.)
        rDstar_lab.SetXYZT(entry.r2_x, entry.r2_y, entry.r2_z, 0.)

        ## 4-momenta        
        pD_lab      = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(entry.q1_x, entry.q1_y, entry.q1_z, mD)
        pDstar_lab  = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(entry.q2_x, entry.q2_y, entry.q2_z, mDstar)

        ## Missing particle momentum and energy
        pDDstar_lab = (pD_lab + pDstar_lab)
        eMiss_lab   = math.sqrt(self.s) - pD_lab.E() - pDstar_lab.E()
        pMiss_lab   = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(- pDDstar_lab.Px(), - pDDstar_lab.Py(), - pDDstar_lab.Pz(), eMiss_lab)

        # Save lab values
        self.r1[:, 0]      = [ rD_lab.T(),     rD_lab.X(),      rD_lab.Y(),      rD_lab.Z()      ]
        self.q1[:, 0]      = [ pD_lab.E(),     pD_lab.Px(),     pD_lab.Py(),     pD_lab.Pz()     ]
        self.r2[:, 0]      = [ rDstar_lab.T(), rDstar_lab.X(),  rDstar_lab.Y(),  rDstar_lab.Z()  ]
        self.q2[:, 0]      = [ pDstar_lab.E(), pDstar_lab.Px(), pDstar_lab.Py(), pDstar_lab.Pz() ]
        self.q3[:, 0]      = [ pMiss_lab.E(),  pMiss_lab.Px(),  pMiss_lab.Py(),  pMiss_lab.Pz()  ]

        # Now we everything to the CM frame
        boost_vector = pDDstar_lab.BoostToCM()
        rD_cm        = ROOT.Math.VectorUtil.boost(rD_lab,      boost_vector)
        pD_cm        = ROOT.Math.VectorUtil.boost(pD_lab,      boost_vector)
        rDstar_cm    = ROOT.Math.VectorUtil.boost(rDstar_lab,  boost_vector)
        pDstar_cm    = ROOT.Math.VectorUtil.boost(pDstar_lab,  boost_vector)
        pDDstar_cm   = ROOT.Math.VectorUtil.boost(pDDstar_lab, boost_vector)
        pMiss_cm     = ROOT.Math.VectorUtil.boost(pMiss_lab,   boost_vector)

        # Save CM values
        self.cm_energy[:,0]    = [ pD_cm.E(),      pDstar_cm.E() ]
        self.cm_masses[:,0]    = [ pDDstar_cm.E(), pMiss_cm.E()  ]
        self.qcm[:,0]          = [ pD_cm.Px(),     pD_cm.Py(),      pD_cm.Pz() ]
        self.rcm1[:,0]         = [ rD_cm.T(),      rD_cm.X(),       rD_cm.Y(),      rD_cm.Z()     ]
        self.rcm2[:,0]         = [ rDstar_cm.T(),  rDstar_cm.X(),   rDstar_cm.Y(),  rDstar_cm.Z() ]

        return