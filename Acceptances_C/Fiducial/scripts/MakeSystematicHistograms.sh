#!/bin/bash

#HIST_DIR=/Users/Chris/CMS/WGamGam/Acceptances/Fiducial/test/Histograms/NoEEEE_NLO_LepGammaGammaFinalElandMuUnblindAll_2015_08_01_ScaleFactors_PDFReweights/

echo python scripts/MakeScaleFactorSystematicHistograms.py 
python scripts/MakeScaleFactorSystematicHistograms.py 


echo python scripts/MakePileUpSystematicHistograms.py 
python scripts/MakePileUpSystematicHistograms.py 


echo python scripts/MakePDFSystematicHistograms.py 
python scripts/MakePDFSystematicHistograms.py 


echo python scripts/MakeFactorizationRenormalizationSystematicHistograms.py 
python scripts/MakeFactorizationRenormalizationSystematicHistograms.py 




