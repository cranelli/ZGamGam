#Python Code for weighting and summing  all
# the Histograms from two root files.
# Example:
# python weightAndAddHistograms.py isr.root fsr.root sum.root

import sys

from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F


# Calculation of ISR and FSR weights
lumi = 19400

# Next to Leading Order
crossSection_isr=0.1592
crossSection_fsr=1.137
#Sum of NLO Weights (Postive and Negative)
# (PtG 500MeV)
N_isr = 154921.6
N_fsr = 1141000.0
#N_isr = 157694.6
#N_fsr = 1123885.0
#N_isr = 992414
#N_fsr = 985000

# Leading Order
#crossSection_isr=0.319
#crossSection_fsr=1.84
#N_isr = 1000310
#N_fsr = 1003920

weight1=crossSection_isr*lumi/N_isr
weight2=crossSection_fsr*lumi/N_fsr


#doCalcC=True

def weightAndAddHistograms(file1Loc, file2Loc, outFileLoc):
    print "weight ISR: ", weight1
    print "weight FSR: ", weight2

    #Iterate over all the Histograms in the first File
    file1 = TFile(file1Loc, 'READ')
    file2 = TFile(file2Loc, 'READ')
    outFile = TFile(outFileLoc, "RECREATE")
    
    list = file1.GetListOfKeys()
    for key in list:
        histName = key.GetName()
        print histName
        hist1 = file1.Get(histName)
        hist2 = file2.Get(histName)
        #Make Sure A Second Histogram of the Same Name Existed
        #print hist2.Class()
        #if hist2.Class() == "TObject" : continue
        weightedSumHist = hist1.Clone(histName)
        weightedSumHist.Reset("CE")

        weightedSumHist.Add(hist1, hist2, weight1, weight2)
        weightedSumHist.Write()
           
if __name__=="__main__":
    weightAndAddHistograms(sys.argv[1], sys.argv[2], sys.argv[3])
    
