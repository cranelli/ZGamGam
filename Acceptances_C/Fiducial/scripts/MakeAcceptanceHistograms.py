#Python Code for looking at the Histograms in a root
#file, and printing the content of each bin.
# Example call:
# python MakeAcceptanceHistograms.py

import sys
import Config

from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
from collections import namedtuple


#Numerator is NO EEEE
#Denominator is NLO Fiducial

##########################################################
#Global Variables (Set from Config File)


RECO_FILE = TFile(Config.RECO_FILE_LOC, 'READ')
GEN_FILE = TFile(Config.GEN_FILE_LOC, 'READ')
ACCEPTANCE_FILE = TFile(Config.ACCEPTANCE_FILE_LOC, "RECREATE")

# Specify the Channels and Hist Types to calculate Acceptances for

RECO_GEN_PAIRS=Config.RECO_GEN_PAIRS
HIST_TYPES=Config.HIST_TYPES

#For Scale Factor Uncertainties
DO_SCALEFACTORS_UPDN=Config.DO_SCALEFACTOR_REWEIGHTS
SFS = Config.SFS

# For Pileup Systematics
DO_PUWEIGHTS_UPDN=Config.DO_PUWEIGHT_REWEIGHTS
PUWEIGHTS = Config.PUWEIGHTS

#Renormalization Reweighting
DO_RENORMALIZATION_REWEIGHT=Config.DO_RENORMALIZATION_REWEIGHT
RENORMALIZATION_VARIATIONS=Config.RENORMALIZATION_VARIATIONS

#Factroization Reweighting
DO_FACTORIZATION_REWEIGHT=Config.DO_FACTORIZATION_REWEIGHT
FACTORIZATION_VARIATIONS=Config.FACTORIZATION_VARIATIONS

# PDF Reweighting
DO_PDFSET_REWEIGHT=Config.DO_PDFSET_REWEIGHT
PDFSET_NAMES=Config.PDFSET_NAMES

DO_PDFEIGENVECTOR_REWEIGHT=Config.DO_PDFEIGENVECTOR_REWEIGHT
PDFEIGENVECTOR_NAME=Config.PDFEIGENVECTOR_NAME
NUM_EIGENVECTORS=Config.NUM_EIGENVECTORS

####################################################################


def MakeAcceptanceHistograms():

    for pair in RECO_GEN_PAIRS:
        for hist_type in HIST_TYPES:
            MakeBasicAcceptanceHistograms(pair, hist_type)
            
            # These Additional Acceptance Histograms,
            # are for use in calculating systematic uncertainties.
            if DO_SCALEFACTORS_UPDN:
                MakeScaleFactorUpDnAcceptanceHistograms(pair, hist_type)
            if DO_PUWEIGHTS_UPDN:
                MakePileUpWeightUpDnAcceptanceHistograms(pair, hist_type)
            if DO_RENORMALIZATION_REWEIGHT:
                MakeRenoramalizationReweightAcceptanceHistograms(pair, hist_type)
            if DO_FACTORIZATION_REWEIGHT:
                MakeFactorizationReweightAcceptanceHistograms(pair, hist_type)
            if DO_PDFSET_REWEIGHT:
                MakePDFSetReweightAcceptanceHistograms(pair, hist_type)
            if DO_PDFEIGENVECTOR_REWEIGHT:
                MakePDFEigenvectorReweightAcceptanceHistograms(pair, hist_type)



# Makes Histograms: Acceptance with Tau as Background, Acceptance with Tau as Signal
def MakeHistograms(recoHistName, genHistName, acceptances_prefix, hist_type):
    print recoHistName
    print genHistName
    recoHist = RECO_FILE.Get(recoHistName)
    genHist = GEN_FILE.Get(genHistName)
    genPlusTauHist = GetGenPlusTauHist(genHistName)
    
    acceptanceTauBkgdName = acceptances_prefix+"_TauBkgd_"+hist_type
    acceptanceTauBkgd = DivideHistograms(recoHist, genHist, acceptanceTauBkgdName)
    acceptanceTauBkgd.Print()
    acceptanceTauBkgd.Write()
    
    acceptanceTauSigName = acceptances_prefix+"_TauSig_"+hist_type
    acceptanceTauSig = DivideHistograms(recoHist, genPlusTauHist, acceptanceTauSigName)
    acceptanceTauSig.Print()
    acceptanceTauSig.Write()
    
    # oneMinusFtau= genHist.GetName()+"_OneMinusFtau"
    # oneMinusFtau = DivideHistograms(genHist, genHistPlustau)
    # oneMinusFtau.Write()


# These are the basic acceptance histograms, gen events are unweighted and
# reco events have the ScaleFactor weights applied.
def MakeBasicAcceptanceHistograms(pair, hist_type):
    #Reco
    recoHistName=pair.reco_channel+"_ScaleFactors_"+hist_type
    #Gen
    genHistName = pair.gen_decay +"_"+hist_type
    #Acceptance Prefix
    acceptances_prefix=pair.reco_channel+"_Acceptance"
    #acceptances_prefix=pair.reco_channel+"_over_"+pair.gen_decay+"_Acceptance"
    MakeHistograms(recoHistName, genHistName, acceptances_prefix, hist_type)


# Calculates the Acceptances for the different possilbe scale factors choices.
# Scale Factor change only effects the Reco Histograms.
# Used in calculating the scale factor systematics.
def MakeScaleFactorUpDnAcceptanceHistograms(pair, hist_type):
    directions = ['UP', 'DN']
    for scalefactor_name in SFS[pair.reco_channel]:
        for dir in directions:
            scalefactor_updn_suffix=scalefactor_name+dir
            #Reco
            recoHistName=pair.reco_channel+"_"+scalefactor_updn_suffix+"_"+hist_type
            #Gen
            genHistName = pair.gen_decay +"_"+hist_type
            #Acceptance Prefix
            acceptance_prefix = pair.reco_channel+"_"+scalefactor_updn_suffix+"_Acceptance"
            MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)

# Calculates the Acceptances for the different possilbe pile up weight choices.
# Used in calculating the pile up systematics.
def MakePileUpWeightUpDnAcceptanceHistograms(pair, hist_type):
    directions = ['UP', 'DN']
    for pu_num in PUWEIGHTS:
        for dir in directions:
            pileup_reweight_suffix="PUWeight"+dir+pu_num
            #Reco
            recoHistName=pair.reco_channel+"_"+pileup_reweight_suffix+"_"+hist_type
            #Gen
            genHistName = pair.gen_decay +"_"+pileup_reweight_suffix+"_"+hist_type
            #Acceptance Prefix
            acceptance_prefix = pair.reco_channel+"_"+pileup_reweight_suffix+"_Acceptance"
            MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)


# Calculates the Acceptances for the different possilbe Renormalization Variations.
# Used in calculating the renormalization systematics.
def MakeRenoramalizationReweightAcceptanceHistograms(pair, hist_type):
    for renormalization_variation in RENORMALIZATION_VARIATIONS:
        renormalization_suffix="Renormalization_" + renormalization_variation
        #Reco
        recoHistName=pair.reco_channel+"_" + renormalization_suffix+"_"+hist_type
        #Gen
        genHistName = pair.gen_decay + "_" + renormalization_suffix+"_"+hist_type
        #Acceptance Prefix
        acceptance_prefix = pair.reco_channel+"_"+renormalization_suffix+"_Acceptance"
        MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)

# Calculates the Acceptances for the different possilbe Factorization Variations.
# Used in calculating the Factorization systematics.
def MakeFactorizationReweightAcceptanceHistograms(pair, hist_type):
    for factorization_variation in FACTORIZATION_VARIATIONS:
        factorization_suffix = "Factorization_"+ factorization_variation
        #Reco
        recoHistName=pair.reco_channel+"_"+ factorization_suffix+"_"+hist_type
        #Gen
        genHistName = pair.gen_decay +"_"+ factorization_suffix + "_" + hist_type
        #Acceptance Prefix
        acceptance_prefix = pair.reco_channel+"_"+factorization_suffix+"_Acceptance"
        MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)

# Calculate the Acceptances for events that have been reweighted to the central value
# of a different PDF set.
def  MakePDFSetReweightAcceptanceHistograms(pair, hist_type):
    for pdf_name in PDFSET_NAMES:
        #pdfset_reweight_suffix=pdf_name+"_PDFReweight"
        pdfset_reweight_suffix = pdf_name
        # Reco
        recoHistName = pair.reco_channel+"_" + pdfset_reweight_suffix + "_"+hist_type
        # Gen
        genHistName = pair.gen_decay+"_"+pdfset_reweight_suffix+"_"+hist_type
        # Acceptance Prefix
        acceptance_prefix = pair.reco_channel+"_"+pdfset_reweight_suffix+"_Acceptance"
        MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)


# Calculate the Acceptances for events that have been reweighted from a PDF set's central value,
# to the same set's eigenvector deviation.
def MakePDFEigenvectorReweightAcceptanceHistograms(pair, hist_type):
    for eigenvector_index in range (0, NUM_EIGENVECTORS):
        pdf_eigenvector_reweight_suffix = PDFEIGENVECTOR_NAME+"_"+str(eigenvector_index)
        # Reco
        recoHistName = pair.reco_channel+"_"+pdf_eigenvector_reweight_suffix+"_"+hist_type
        # Gen
        genHistName = pair.gen_decay+"_"+pdf_eigenvector_reweight_suffix+"_"+hist_type
        # Acceptance Prefix
        acceptance_prefix = pair.reco_channel+"_"+pdf_eigenvector_reweight_suffix+"_Acceptance"
        MakeHistograms(recoHistName, genHistName, acceptance_prefix, hist_type)

#From the Generator Name, also gets the Generator Histogram for the same decay through the tau.
# Ie Muon and TauToMuon decay.  Then returns the sume of the two histograms.
def GetGenPlusTauHist(genHistName):
    print genHistName
    genHist = GEN_FILE.Get(genHistName)
    
    
    genTauHistName = "TauTo" + genHistName
    genTauHist = GEN_FILE.Get(genTauHistName)
    
    
    genPlusTauHistName = genHistName+"_PlusTau"
    genPlusTauHist = genHist.Clone(genPlusTauHistName)
    genPlusTauHist.Add(genTauHist)

    return genPlusTauHist


# Given Two Histogrmas, divide them, and return the new histogram with the given name.
def DivideHistograms(num_hist, den_hist, divide_hist_name):
    divideHist = num_hist.Clone(divide_hist_name)
    divideHist.Divide(den_hist)
    return divideHist





if __name__=="__main__":
    MakeAcceptanceHistograms()
