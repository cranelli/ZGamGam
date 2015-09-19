#!/bin/bash

#
# Script for running over Final Skim RECO ntuples, 
# to create Histograms used for Acceptance and
# Systematic Calculations.
# Example Call:
# ./scripts/make_reco_hists.sh
#

RECO_NTUPLES_DIR=/data/users/cranelli/WGamGam/ReFilterFinalNtuple/NLO_LepGammaGammaFinalElandMuUnblindAllPtG500MeV_2015_09_16_ScaleFactors_PDFReweights/
HISTS_OUTDIR=test/Histograms/NoEEEE_NLO_LepGammaGammaFinalElandMuUnblindAllPtG500MeV_2015_09_16_ScaleFactors_PDFReweights/

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

# Add the ISR and FSR histograms together, for a Weighted Total Histogram
echo python scripts/weightAndAddHistograms.py $HISTS_OUTDIR"/"job_NLO_WAA_ISR_RecoHistograms.root $HISTS_OUTDIR"/"job_NLO_WAA_FSR_RecoHistograms.root $HISTS_OUTDIR"/"job_NLO_WAA_WeightedTotal_RecoHistograms.root
python scripts/weightAndAddHistograms.py $HISTS_OUTDIR"/"job_NLO_WAA_ISR_RecoHistograms.root $HISTS_OUTDIR"/"job_NLO_WAA_FSR_RecoHistograms.root $HISTS_OUTDIR"/"job_NLO_WAA_WeightedTotal_RecoHistograms.root