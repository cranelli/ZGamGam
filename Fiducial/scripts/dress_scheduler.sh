#!/bin/bash

#
# Script for running over ntuples, to create 
# common fiducial skim.
# Example Call:
# ./scripts/dress_scheduler.sh
#

DRESS_NTUPLES_DIR=/data/users/cranelli/WGamGam/NLO_ggNtuples/CustomForDress/TruthSkim/
DRESS_OUTDIR=/data/users/cranelli/WGamGam/Acceptances/Dressed/CommonFiducial_NLO_wMT_Dress500MeV_Skim/

#CommonFiducial/Dress/CommonFiducial_NLO_wMT_Dress500MeV_Skim

# Run makefile
make dress

#NLO Samples
for prefix in ISR FSR
 do
    in_file="job_NLO_WAA_"$prefix"_PtG500MeV_TruthSkim.root" 
    echo "./test/CommonFiducialDressLeptonSkim.exe" $DRESS_NTUPLES_DIR"/"$in_file $DRESS_OUTDIR"/job_NLO_WAA_"$prefix"/Job_0000/tree.root"

    ./test/CommonFiducialDressLeptonSkim.exe $DRESS_NTUPLES_DIR"/"$in_file $DRESS_OUTDIR"/job_NLO_WAA_"$prefix"/Job_0000/tree.root"
done
