#Python Code for weighting and summing  all
# the Histograms from two root files.
# Example:
# python scaleHistograms.py isr.root fsr.root sum.root

import sys

from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F


# Luminosity, NLO cross section, and Number of weighted (NLO) events,
# for the NLO Zgg sample.

LUMI = 19700

CROSS_SECTION = 0.7518
#Sum of NLO Weights (Postive and Negative)
# (PtG 500MeV)
N_WEIGHTED = 375828

weight=CROSS_SECTION*LUMI/N_WEIGHTED



def scaleHistograms(file1Loc, outFileLoc):
    print "weight: ", weight

    #Iterate over all the Histograms in the first File
    file1 = TFile(file1Loc, 'READ')
    #file2 = TFile(file2Loc, 'READ')
    outFile = TFile(outFileLoc, "RECREATE")
    
    dir = file1.Get("ggNtuplizer")
    list = dir.GetListOfKeys()
    for key in list:
        histName = key.GetName()
        #print histName
        hist1 = dir.Get(histName)
        #print hist1.ClassName()
        
        if hist1.ClassName() == "TTree" : continue
        
        #print hist1.Class()
        hist1.Scale(weight)
        hist1.Write()
           
if __name__=="__main__":
    scaleHistograms(sys.argv[1], sys.argv[3])
    
