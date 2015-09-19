#
# Configuration File, to store all the information
# needed to run the Acceptance and Systematic scripts
#

from collections import namedtuple

RECO_FILE_LOC = "test/Histograms/NoEEEE_NLO_LepGammaGammaFinalElandMuUnblindAllPtG500MeV_2015_09_16_ScaleFactors_PDFReweights/job_NLO_WAA_WeightedTotal_RecoHistograms.root"

GEN_DIR = "test/Histograms/CommonFiducial_NLO_wMT_Dress500MeV_Skim_PUWeights_PDFReweights/"
GEN_FILE_LOC = GEN_DIR + "job_NLO_WAA_WeightedTotal_GenHistograms.root"

HIST_DIR =  "test/Histograms/NoEEEE_NLO_LepGammaGammaFinalElandMuUnblindAll_2015_09_14_ScaleFactors_PDFReweights/"

ACCEPTANCE_FILE_LOC= HIST_DIR + "Acceptances.root"

PILEUP_SYSTEMATICS_FILE_LOC = HIST_DIR + "PileUpSystematics.root"

SCALEFACTOR_SYSTEMATICS_FILE_LOC = HIST_DIR + "ScaleFactorSystematics.root"

PDF_SYSTEMATICS_FILE_LOC = HIST_DIR + "PDFSystematics.root"

FACTORIZATION_RENORMALIZATION_SYSTEMATICS_FILE_LOC = HIST_DIR + "FactorizationRenormalizationSystematics.root"

TOTAL_SYSTEMATICS_FILE_LOC= HIST_DIR + "TotalSystematics.root"

# Specify the Channels and Hist Types to calculate Acceptances for  
CHANNELS=["MuonChannel", "ElectronChannel"]

ChannelPair = namedtuple('ChannelPair','reco_channel gen_decay')
RECO_GEN_PAIRS=[ChannelPair('MuonChannel', 'MuonDecay'), ChannelPair('ElectronChannel', 'ElectronDecay')]

HIST_TYPES=["Count", "Pt", "Category_Pt"]

#For Scale Factor Uncertainties 
DO_SCALEFACTOR_REWEIGHTS=True
# Channel name must match choice in RECO_GEN_PAIRS                                                                                           
SFS={ 'MuonChannel': ["mu_trigSF", "mu_isoSF", "mu_idSF", "ph_idSF"],
    'ElectronChannel': ["el_trigSF", "ph_idSF", "ph_evetoSF"] }

# For Pileup Systematics
DO_PUWEIGHT_REWEIGHTS=True
PUWEIGHTS = ["5"]

#Renormalization Reweighting                                                                                                                 
DO_RENORMALIZATION_REWEIGHT=True
RENORMALIZATION_VARIATIONS=["Half", "Double"]

#Factroization Reweighting                                                                                                                   
DO_FACTORIZATION_REWEIGHT=True
FACTORIZATION_VARIATIONS=["Half", "Double"]

# PDF Reweighting                                                                                                                            
DO_PDFSET_REWEIGHT=True
PDFSET_NAMES=['NNPDF30_nlo_nf_5_pdfas', 'CT10nlo', 'MSTW2008nlo68cl']
ORIG_PDF_NAME='NNPDF30_nlo_nf_5_pdfas'
#PDFSET_NAMES =[ "cteq6l1", "MSTW2008lo68cl", "cteq66" ] 

DO_PDFEIGENVECTOR_REWEIGHT=True
PDFEIGENVECTOR_NAME= 'NNPDF30_nlo_nf_5_pdfas'
#PDFEIGENVECTOR_NAME= "cteq66"                                                                                                               
NUM_EIGENVECTORS=103


#Acceptances
ACCEPTANCE_TYPES=["TauSig", "TauBkgd"]
