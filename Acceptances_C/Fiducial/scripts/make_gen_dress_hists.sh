#!/bin/bash

#
# Script for running over Final Skim RECO ntuples, 
# to create Histograms used for Acceptance and
# Systematic Calculations.
# Example Call:
# ./scripts/make_gen_dressed_hists.sh
#

GEN_DRESS_SKIM_DIR=/data/users/cranelli/ZGamGam/Fiducial/Dressed/CommonFiducial_NLO_Skim_2015_09_22_PDFReweights/

HISTS_OUTDIR=test/Histograms/CommonFiducial_NLO_Skim_2015_09_22/

#CommonFiducial/Dress/CommonFiducial_NLO_wMT_Dress500MeV_Skim

# Run makefile
make gen

# Run over all samples in the directory
for sample in $( ls $GEN_DRESS_SKIM_DIR )
do
    in_file=$GEN_DRESS_SKIM_DIR$sample"/tree.root"
    echo ./test/GenHistograms.exe $in_file $HISTS_OUTDIR"/"$sample"_GenHistograms.root"
    ./test/GenHistograms.exe $in_file $HISTS_OUTDIR"/"$sample"_GenHistograms.root"
done
