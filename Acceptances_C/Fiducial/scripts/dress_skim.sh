#!/bin/bash

#
# Script for running over ntuples, to create 
# common fiducial skim.
# Example Call:
# ./scripts/dress_skim.sh
#

NTUPLES_DIR=/data/users/cranelli/ZGamGam/NLO_ggNtuples/
DRESS_OUTDIR=/data/users/cranelli/ZGamGam/Fiducial/Dressed/CommonFiducial_NLO_Skim/

#CommonFiducial/Dress/CommonFiducial_NLO_wMT_Dress500MeV_Skim

# Run makefile
make dress

#NLO Samples

in_file=llaa_nlo_ggNtuple_part1.root
echo "./test/CommonFiducialDressLeptonSkim.exe" $NTUPLES_DIR"/"$in_file $DRESS_OUTDIR"/job_llaa_nlo/Job_0000/tree.root"

./test/CommonFiducialDressLeptonSkim.exe" $NTUPLES_DIR"/"$in_file $DRESS_OUTDIR"/job_llaa_nlo/Job_0000/tree.root
