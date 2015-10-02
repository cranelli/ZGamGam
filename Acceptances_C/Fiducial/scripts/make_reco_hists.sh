#!/bin/bash

#
# Script for running over Final Skim RECO ntuples, 
# to create Histograms used for Acceptance and
# Systematic Calculations.
# Example Call:
# ./scripts/make_reco_hists.sh
#

RECO_NTUPLES_DIR=/data/users/cranelli/ZGamGam/ReFilterFinalNtuple/NLO_LepLepGammaGammaFinalElElandMuMuUnblindAll_2015_10_01_ScaleFactors_PDFReweights/
HISTS_OUTDIR=test/Histograms/No_EEEE_NLO_LepLepGammaGammaFinalElElandMuMuUnblindAll_2015_10_01_ScaleFactors_PDFReweights/

#CommonFiducial/Dress/CommonFiducial_NLO_wMT_Dress500MeV_Skim

# Run makefile
make reco

# Run over all samples in the directory
for sample in $( ls $RECO_NTUPLES_DIR )
do
    in_file=$RECO_NTUPLES_DIR$sample"/tree.root"
    echo ./test/RecoHistograms.exe $in_file $HISTS_OUTDIR"/"$sample"_RecoHistograms.root"
    ./test/RecoHistograms.exe $in_file $HISTS_OUTDIR"/"$sample"_RecoHistograms.root"
    
done