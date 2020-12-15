# Here I put all the different histograms I want to make
#
# Author:       Daniel Winney (2020)
# Affiliation:  Joint Physics Analysis Center (JPAC)
# Email:        dwinney@iu.edu
# ---------------------------------------------------------------------------

try:
    import ROOT
except:
    exit()

def make_histograms(file):
    h1 = ROOT.TH2F('rel_mom_v_sep', 'Lab Frame', 10, 0, 50, 10, 0, 50)

    for entry in file.lab:  
        rel_mom_vs_sep(h1, entry)

    # h1.Draw("COLZ")
    h1.Write()
    
def rel_mom_vs_sep(h, entry):
    # calculate relative seperation
    xD      = ROOT.TVector3(entry.D_x, entry.D_y, entry.D_z)
    xDstar  = ROOT.TVector3(entry.Dstar_x, entry.Dstar_y, entry.Dstar_z)
    xRel = xD - xDstar
    relSep = xRel.Mag()
    relSep *= 0.1975 # convert to fm 

    # relative momenta
    pD      = ROOT.TVector3(entry.D_px, entry.D_py, entry.D_pz)
    pDstar  = ROOT.TVector3(entry.Dstar_px, entry.Dstar_py, entry.Dstar_pz)
    pRel = pD - pDstar
    relMom = pRel.Mag() 

    # plot style settings
    h.GetXaxis().SetTitle("Rel. Seperation [fm]")
    h.GetYaxis().SetTitle("Rel. Momentum [GeV]")
    h.SetOption("COLZ")     # colored bins
    h.SetStats(0)           # remove stats box
    h.Fill(relSep, relMom)