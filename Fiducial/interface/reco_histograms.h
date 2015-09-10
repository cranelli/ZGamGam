//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  9 16:17:52 2015 by ROOT version 5.34/30
// from TTree EventTree/EventTree
// found on file: /data/users/cranelli/WGamGam/ReFilterFinalNtuple/NLO_LepGammaGammaFinalElandMuUnblindAll_2015_08_01_ScaleFactors_PDFReweights/job_NLO_WAA_ISR/tree.root
//////////////////////////////////////////////////////////

#ifndef reco_histograms_h
#define reco_histograms_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "mc_particle_data.h"
#include "histogram_builder.h"
#include "cut_values.h"

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <map>
#include <string>

// Fixed size dimensions of array or collections stored in the TTree if any.

class RecoHistograms {
 private :
  string SelectChannel();
  string SelectPhotonsLocation();
  
 public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain


   //My Stuff
   HistogramBuilder histogram_builder_;
   void MakeHistograms(string prefix, double weight);   
   void MakeBasicHistograms(string channel_type);   
   void MakeUnweightedHistograms(string channel_type);

   // Declaration of leaf types
   vector<float>   *LHEWeight_weights;
   vector<string>  *LHEWeight_ids;
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nVtxBS;
   Float_t         pdf[7];
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcPt;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   vector<int>     *nPU;
   vector<float>   *puTrue;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfType01MET;
   Float_t         pfType01METPhi;
   Float_t         pfType01METsumEt;
   Float_t         pfType01METmEtSig;
   Float_t         pfType01METSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Float_t         rho2012;
   Int_t           el_n;
   Int_t           mu_n;
   Int_t           ph_n;
   Int_t           jet_n;
   Int_t           vtx_n;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_sceta;
   vector<float>   *el_phi;
   vector<float>   *el_e;
   vector<float>   *el_pt_uncorr;
   vector<float>   *el_e_uncorr;
   vector<float>   *el_pfiso30;
   vector<float>   *el_pfiso40;
   vector<int>     *el_charge;
   vector<bool>    *el_triggerMatch;
   vector<bool>    *el_hasMatchedConv;
   vector<bool>    *el_passTight;
   vector<bool>    *el_passMedium;
   vector<bool>    *el_passLoose;
   vector<bool>    *el_passVeryLoose;
   vector<bool>    *el_passTightTrig;
   vector<bool>    *el_passMvaTrig;
   vector<bool>    *el_passMvaNonTrig;
   vector<bool>    *el_passMvaTrigNoIso;
   vector<bool>    *el_passMvaNonTrigNoIso;
   vector<bool>    *el_passMvaTrigOnlyIso;
   vector<bool>    *el_passMvaNonTrigOnlyIso;
   vector<bool>    *el_truthMatch_el;
   vector<float>   *el_truthMinDR_el;
   vector<float>   *el_truthMatchPt_el;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_e;
   vector<float>   *mu_pt_uncorr;
   vector<float>   *mu_eta_uncorr;
   vector<float>   *mu_phi_uncorr;
   vector<float>   *mu_e_uncorr;
   vector<float>   *mu_corrIso;
   vector<int>     *mu_charge;
   vector<bool>    *mu_triggerMatch;
   vector<bool>    *mu_triggerMatchDiMu;
   vector<bool>    *mu_passTight;
   vector<bool>    *mu_passTightNoIso;
   vector<bool>    *mu_passTightNoD0;
   vector<bool>    *mu_passTightNoIsoNoD0;
   vector<bool>    *mu_truthMatch;
   vector<float>   *mu_truthMinDR;
   vector<float>   *ph_pt;
   vector<float>   *ph_eta;
   vector<float>   *ph_sceta;
   vector<float>   *ph_phi;
   vector<float>   *ph_e;
   vector<float>   *ph_scE;
   vector<float>   *ph_pt_uncorr;
   vector<float>   *ph_HoverE;
   vector<float>   *ph_HoverE12;
   vector<float>   *ph_sigmaIEIE;
   vector<float>   *ph_chIsoCorr;
   vector<float>   *ph_neuIsoCorr;
   vector<float>   *ph_phoIsoCorr;
   vector<bool>    *ph_eleVeto;
   vector<bool>    *ph_hasPixSeed;
   vector<bool>    *ph_isConv;
   vector<bool>    *ph_passTight;
   vector<bool>    *ph_passMedium;
   vector<bool>    *ph_passLoose;
   vector<bool>    *ph_passLooseNoSIEIE;
   vector<bool>    *ph_passHOverELoose;
   vector<bool>    *ph_passHOverEMedium;
   vector<bool>    *ph_passHOverETight;
   vector<bool>    *ph_passSIEIELoose;
   vector<bool>    *ph_passSIEIEMedium;
   vector<bool>    *ph_passSIEIETight;
   vector<bool>    *ph_passChIsoCorrLoose;
   vector<bool>    *ph_passChIsoCorrMedium;
   vector<bool>    *ph_passChIsoCorrTight;
   vector<bool>    *ph_passNeuIsoCorrLoose;
   vector<bool>    *ph_passNeuIsoCorrMedium;
   vector<bool>    *ph_passNeuIsoCorrTight;
   vector<bool>    *ph_passPhoIsoCorrLoose;
   vector<bool>    *ph_passPhoIsoCorrMedium;
   vector<bool>    *ph_passPhoIsoCorrTight;
   vector<bool>    *ph_truthMatch_el;
   vector<float>   *ph_truthMinDR_el;
   vector<float>   *ph_truthMatchPt_el;
   vector<bool>    *ph_truthMatch_ph;
   vector<float>   *ph_truthMinDR_ph;
   vector<float>   *ph_truthMatchPt_ph;
   vector<int>     *ph_truthMatchMotherPID_ph;
   vector<bool>    *ph_IsEB;
   vector<bool>    *ph_IsEE;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_JECUnc;
   vector<float>   *jet_e;
   vector<bool>    *jet_PUIDLoose;
   vector<bool>    *jet_PUIDMedium;
   vector<bool>    *jet_PUIDTight;
   vector<float>   *jet_CSV;
   vector<int>     *jet_genIndex;
   vector<float>   *jet_genPt;
   vector<float>   *jet_genEta;
   vector<float>   *jet_genPhi;
   vector<float>   *jet_genE;
   Float_t         PUWeight;
   Float_t         PUWeightDN5;
   Float_t         PUWeightDN10;
   Float_t         PUWeightUP5;
   Float_t         PUWeightUP6;
   Float_t         PUWeightUP7;
   Float_t         PUWeightUP8;
   Float_t         PUWeightUP9;
   Float_t         PUWeightUP10;
   Float_t         PUWeightUP11;
   Float_t         PUWeightUP12;
   Float_t         PUWeightUP13;
   Float_t         PUWeightUP14;
   Float_t         PUWeightUP15;
   Float_t         PUWeightUP16;
   Float_t         PUWeightUP17;
   Bool_t          passTrig_ele27WP80;
   Bool_t          passTrig_mu24eta2p1;
   Bool_t          passTrig_mu24;
   Bool_t          passTrig_mu17_mu8;
   Bool_t          passTrig_mu17_Tkmu8;
   Bool_t          passTrig_ele17_ele8_9;
   Bool_t          passTrig_ele17_ele8_22;
   Bool_t          isBlinded;
   Float_t         EventWeight;
   vector<bool>    *ph_hasMatchedEle;
   Int_t           mu_pt25_n;
   Int_t           mu_passtrig_n;
   Int_t           mu_passtrig25_n;
   Int_t           el_pt25_n;
   Int_t           el_passtrig_n;
   Int_t           el_passtrig28_n;
   Int_t           ph_mediumNoSIEIE_n;
   Int_t           ph_medium_n;
   Int_t           ph_mediumNoEleVeto_n;
   Int_t           ph_mediumNoSIEIENoEleVeto_n;
   Int_t           ph_mediumNoChIsoNoEleVeto_n;
   Int_t           ph_mediumNoNeuIsoNoEleVeto_n;
   Int_t           ph_mediumNoPhoIsoNoEleVeto_n;
   Int_t           ph_mediumNoIso_n;
   Int_t           ph_mediumNoChIso_n;
   Int_t           ph_mediumNoNeuIso_n;
   Int_t           ph_mediumNoPhoIso_n;
   Int_t           ph_mediumNoChIsoNoNeuIso_n;
   Int_t           ph_mediumNoChIsoNoPhoIso_n;
   Int_t           ph_mediumNoNeuIsoNoPhoIso_n;
   Int_t           ph_iso533_n;
   Int_t           ph_iso855_n;
   Int_t           ph_iso1077_n;
   Int_t           ph_iso1299_n;
   Int_t           ph_iso151111_n;
   Int_t           ph_iso201616_n;
   vector<bool>    *ph_trigMatch_el;
   vector<float>   *ph_elMinDR;
   Float_t         leadPhot_pt;
   Float_t         sublPhot_pt;
   Float_t         leadPhot_lepDR;
   Float_t         sublPhot_lepDR;
   Float_t         ph_phDR;
   Float_t         phPhot_lepDR;
   Float_t         leadPhot_lepDPhi;
   Float_t         sublPhot_lepDPhi;
   Float_t         ph_phDPhi;
   Float_t         phPhot_lepDPhi;
   Float_t         dphi_met_lep1;
   Float_t         dphi_met_lep2;
   Float_t         dphi_met_ph1;
   Float_t         dphi_met_ph2;
   Float_t         mt_lep_met;
   Float_t         mt_lepph1_met;
   Float_t         mt_lepph2_met;
   Float_t         mt_lepphph_met;
   Float_t         m_leplep;
   Float_t         m_mumu;
   Float_t         m_elel;
   Float_t         m_leplep_uncorr;
   Float_t         m_lepph1;
   Float_t         m_lepph2;
   Float_t         m_lep2ph1;
   Float_t         m_lep2ph2;
   Float_t         m_lepphlead;
   Float_t         m_lepphsubl;
   Float_t         m_lep2phlead;
   Float_t         m_lep2phsubl;
   Float_t         m_leplepph;
   Float_t         m_leplepphph;
   Float_t         m_leplepph1;
   Float_t         m_leplepph2;
   Float_t         m_lepphph;
   Float_t         m_leplepZ;
   Float_t         m_3lep;
   Float_t         m_4lep;
   Float_t         pt_leplep;
   Float_t         pt_lepph1;
   Float_t         pt_lepph2;
   Float_t         pt_lepphph;
   Float_t         pt_leplepph;
   Float_t         pt_secondLepton;
   Float_t         pt_thirdLepton;
   Float_t         leadPhot_leadLepDR;
   Float_t         leadPhot_sublLepDR;
   Float_t         sublPhot_leadLepDR;
   Float_t         sublPhot_sublLepDR;
   Float_t         dr_ph1_leadLep;
   Float_t         dr_ph1_sublLep;
   Float_t         dr_ph2_leadLep;
   Float_t         dr_ph2_sublLep;
   Float_t         dphi_ph1_leadLep;
   Float_t         dphi_ph1_sublLep;
   Float_t         dphi_ph2_leadLep;
   Float_t         dphi_ph2_sublLep;
   Float_t         m_ph1_ph2;
   Float_t         dr_ph1_ph2;
   Float_t         dphi_ph1_ph2;
   Float_t         pt_ph1_ph2;
   Float_t         m_leadLep_ph1_ph2;
   Float_t         m_leadLep_ph1;
   Float_t         m_leadLep_ph2;
   Float_t         pt_leadph12;
   Float_t         pt_sublph12;
   Float_t         eta_leadph12;
   Float_t         eta_sublph12;
   Float_t         hasPixSeed_leadph12;
   Float_t         hasPixSeed_sublph12;
   Float_t         sieie_leadph12;
   Float_t         sieie_sublph12;
   Float_t         chIsoCorr_leadph12;
   Float_t         chIsoCorr_sublph12;
   Float_t         neuIsoCorr_leadph12;
   Float_t         neuIsoCorr_sublph12;
   Float_t         phoIsoCorr_leadph12;
   Float_t         phoIsoCorr_sublph12;
   Bool_t          isEB_leadph12;
   Bool_t          isEB_sublph12;
   Bool_t          isEE_leadph12;
   Bool_t          isEE_sublph12;
   Bool_t          truthMatchPh_leadph12;
   Bool_t          truthMatchPh_sublph12;
   Bool_t          truthMatchPhMomPID_leadph12;
   Bool_t          truthMatchPhMomPID_sublph12;
   Float_t         m_nearestToZ;
   Float_t         m_minZdifflepph;
   Int_t           truelep_n;
   Int_t           trueph_n;
   Int_t           trueph_wmother_n;
   Int_t           truegenph_n;
   Int_t           truegenphpt15_n;
   Int_t           truegenphpt15WZMom;
   Int_t           truegenphpt15LepMom_n;
   Int_t           truegenphpt15QMom_n;
   vector<float>   *truelep_pt;
   vector<float>   *truelep_eta;
   vector<float>   *truelep_phi;
   vector<float>   *truelep_e;
   vector<bool>    *truelep_isElec;
   vector<bool>    *truelep_isMuon;
   vector<int>     *truelep_motherPID;
   vector<float>   *trueph_pt;
   vector<float>   *trueph_eta;
   vector<float>   *trueph_phi;
   vector<int>     *trueph_motherPID;
   vector<int>     *trueph_parentage;
   vector<float>   *trueph_nearestLepDR;
   vector<float>   *trueph_nearestQrkDR;
   vector<float>   *trueW_pt;
   vector<float>   *trueW_eta;
   vector<float>   *trueW_phi;
   vector<float>   *trueW_e;
   Float_t         trueleadlep_pt;
   Float_t         truesubllep_pt;
   Float_t         true_m_leplep;
   Float_t         trueleadlep_leadPhotDR;
   Float_t         trueleadlep_sublPhotDR;
   Float_t         truesubllep_leadPhotDR;
   Float_t         truesubllep_sublPhotDR;
   Float_t         truephph_dr;
   Float_t         truephph_dphi;
   Float_t         truephph_m;
   Float_t         truelepphph_m;
   Float_t         el_trigSF;
   Float_t         el_trigSFUP;
   Float_t         el_trigSFDN;
   Float_t         ph_idSF;
   Float_t         ph_idSFUP;
   Float_t         ph_idSFDN;
   Float_t         ph_evetoSF;
   Float_t         ph_evetoSFUP;
   Float_t         ph_evetoSFDN;
   Float_t         mu_trigSF;
   Float_t         mu_trigSFUP;
   Float_t         mu_trigSFDN;
   Float_t         mu_isoSF;
   Float_t         mu_isoSFUP;
   Float_t         mu_isoSFDN;
   Float_t         mu_idSF;
   Float_t         mu_idSFUP;
   Float_t         mu_idSFDN;
   vector<double>  *xfx_first_NNPDF30_nlo_nf_5_pdfas;
   vector<double>  *xfx_second_NNPDF30_nlo_nf_5_pdfas;
   vector<double>  *xfx_first_CT10nlo;
   vector<double>  *xfx_second_CT10nlo;
   vector<double>  *xfx_first_MSTW2008nlo68cl;
   vector<double>  *xfx_second_MSTW2008nlo68cl;

   //My leaf types

   //PileUp Reweights                                                    
   map<string, Float_t>  PU_Reweights;
   // PDF Reweights - Contains two parton distribution functions, one for each proton                          
   map<string, pair<vector<double> *, vector<double> * > > PDF_Reweights;   

   // List of branches
   TBranch        *b_LHEWeight_weights;   //!
   TBranch        *b_LHEWeight_ids;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfType01MET;   //!
   TBranch        *b_pfType01METPhi;   //!
   TBranch        *b_pfType01METsumEt;   //!
   TBranch        *b_pfType01METmEtSig;   //!
   TBranch        *b_pfType01METSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_ph_n;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_vtx_n;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_sceta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_pt_uncorr;   //!
   TBranch        *b_el_e_uncorr;   //!
   TBranch        *b_el_pfiso30;   //!
   TBranch        *b_el_pfiso40;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_triggerMatch;   //!
   TBranch        *b_el_hasMatchedConv;   //!
   TBranch        *b_el_passTight;   //!
   TBranch        *b_el_passMedium;   //!
   TBranch        *b_el_passLoose;   //!
   TBranch        *b_el_passVeryLoose;   //!
   TBranch        *b_el_passTightTrig;   //!
   TBranch        *b_el_passMvaTrig;   //!
   TBranch        *b_el_passMvaNonTrig;   //!
   TBranch        *b_el_passMvaTrigNoIso;   //!
   TBranch        *b_el_passMvaNonTrigNoIso;   //!
   TBranch        *b_el_passMvaTrigOnlyIso;   //!
   TBranch        *b_el_passMvaNonTrigOnlyIso;   //!
   TBranch        *b_el_truthMatch_el;   //!
   TBranch        *b_el_truthMinDR_el;   //!
   TBranch        *b_el_truthMatchPt_el;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_pt_uncorr;   //!
   TBranch        *b_mu_eta_uncorr;   //!
   TBranch        *b_mu_phi_uncorr;   //!
   TBranch        *b_mu_e_uncorr;   //!
   TBranch        *b_mu_corrIso;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_triggerMatch;   //!
   TBranch        *b_mu_triggerMatchDiMu;   //!
   TBranch        *b_mu_passTight;   //!
   TBranch        *b_mu_passTightNoIso;   //!
   TBranch        *b_mu_passTightNoD0;   //!
   TBranch        *b_mu_passTightNoIsoNoD0;   //!
   TBranch        *b_mu_truthMatch;   //!
   TBranch        *b_mu_truthMinDR;   //!
   TBranch        *b_ph_pt;   //!
   TBranch        *b_ph_eta;   //!
   TBranch        *b_ph_sceta;   //!
   TBranch        *b_ph_phi;   //!
   TBranch        *b_ph_e;   //!
   TBranch        *b_ph_scE;   //!
   TBranch        *b_ph_pt_uncorr;   //!
   TBranch        *b_ph_HoverE;   //!
   TBranch        *b_ph_HoverE12;   //!
   TBranch        *b_ph_sigmaIEIE;   //!
   TBranch        *b_ph_chIsoCorr;   //!
   TBranch        *b_ph_neuIsoCorr;   //!
   TBranch        *b_ph_phoIsoCorr;   //!
   TBranch        *b_ph_eleVeto;   //!
   TBranch        *b_ph_hasPixSeed;   //!
   TBranch        *b_ph_isConv;   //!
   TBranch        *b_ph_passTight;   //!
   TBranch        *b_ph_passMedium;   //!
   TBranch        *b_ph_passLoose;   //!
   TBranch        *b_ph_passLooseNoSIEIE;   //!
   TBranch        *b_ph_passHOverELoose;   //!
   TBranch        *b_ph_passHOverEMedium;   //!
   TBranch        *b_ph_passHOverETight;   //!
   TBranch        *b_ph_passSIEIELoose;   //!
   TBranch        *b_ph_passSIEIEMedium;   //!
   TBranch        *b_ph_passSIEIETight;   //!
   TBranch        *b_ph_passChIsoCorrLoose;   //!
   TBranch        *b_ph_passChIsoCorrMedium;   //!
   TBranch        *b_ph_passChIsoCorrTight;   //!
   TBranch        *b_ph_passNeuIsoCorrLoose;   //!
   TBranch        *b_ph_passNeuIsoCorrMedium;   //!
   TBranch        *b_ph_passNeuIsoCorrTight;   //!
   TBranch        *b_ph_passPhoIsoCorrLoose;   //!
   TBranch        *b_ph_passPhoIsoCorrMedium;   //!
   TBranch        *b_ph_passPhoIsoCorrTight;   //!
   TBranch        *b_ph_truthMatch_el;   //!
   TBranch        *b_ph_truthMinDR_el;   //!
   TBranch        *b_ph_truthMatchPt_el;   //!
   TBranch        *b_ph_truthMatch_ph;   //!
   TBranch        *b_ph_truthMinDR_ph;   //!
   TBranch        *b_ph_truthMatchPt_ph;   //!
   TBranch        *b_ph_truthMatchMotherPID_ph;   //!
   TBranch        *b_ph_IsEB;   //!
   TBranch        *b_ph_IsEE;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_JECUnc;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_PUIDLoose;   //!
   TBranch        *b_jet_PUIDMedium;   //!
   TBranch        *b_jet_PUIDTight;   //!
   TBranch        *b_jet_CSV;   //!
   TBranch        *b_jet_genIndex;   //!
   TBranch        *b_jet_genPt;   //!
   TBranch        *b_jet_genEta;   //!
   TBranch        *b_jet_genPhi;   //!
   TBranch        *b_jet_genE;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_PUWeightDN5;   //!
   TBranch        *b_PUWeightDN10;   //!
   TBranch        *b_PUWeightUP5;   //!
   TBranch        *b_PUWeightUP6;   //!
   TBranch        *b_PUWeightUP7;   //!
   TBranch        *b_PUWeightUP8;   //!
   TBranch        *b_PUWeightUP9;   //!
   TBranch        *b_PUWeightUP10;   //!
   TBranch        *b_PUWeightUP11;   //!
   TBranch        *b_PUWeightUP12;   //!
   TBranch        *b_PUWeightUP13;   //!
   TBranch        *b_PUWeightUP14;   //!
   TBranch        *b_PUWeightUP15;   //!
   TBranch        *b_PUWeightUP16;   //!
   TBranch        *b_PUWeightUP17;   //!
   TBranch        *b_passTrig_ele27WP80;   //!
   TBranch        *b_passTrig_mu24eta2p1;   //!
   TBranch        *b_passTrig_mu24;   //!
   TBranch        *b_passTrig_mu17_mu8;   //!
   TBranch        *b_passTrig_mu17_Tkmu8;   //!
   TBranch        *b_passTrig_ele17_ele8_9;   //!
   TBranch        *b_passTrig_ele17_ele8_22;   //!
   TBranch        *b_isBlinded;   //!
   TBranch        *b_EventWeight;   //!
   TBranch        *b_ph_hasMatchedEle;   //!
   TBranch        *b_mu_pt25_n;   //!
   TBranch        *b_mu_passtrig_n;   //!
   TBranch        *b_mu_passtrig25_n;   //!
   TBranch        *b_el_pt25_n;   //!
   TBranch        *b_el_passtrig_n;   //!
   TBranch        *b_el_passtrig28_n;   //!
   TBranch        *b_ph_mediumNoSIEIE_n;   //!
   TBranch        *b_ph_medium_n;   //!
   TBranch        *b_ph_mediumNoEleVeto_n;   //!
   TBranch        *b_ph_mediumNoSIEIENoEleVeto_n;   //!
   TBranch        *b_ph_mediumNoChIsoNoEleVeto_n;   //!
   TBranch        *b_ph_mediumNoNeuIsoNoEleVeto_n;   //!
   TBranch        *b_ph_mediumNoPhoIsoNoEleVeto_n;   //!
   TBranch        *b_ph_mediumNoIso_n;   //!
   TBranch        *b_ph_mediumNoChIso_n;   //!
   TBranch        *b_ph_mediumNoNeuIso_n;   //!
   TBranch        *b_ph_mediumNoPhoIso_n;   //!
   TBranch        *b_ph_mediumNoChIsoNoNeuIso_n;   //!
   TBranch        *b_ph_mediumNoChIsoNoPhoIso_n;   //!
   TBranch        *b_ph_mediumNoNeuIsoNoPhoIso_n;   //!
   TBranch        *b_ph_iso533_n;   //!
   TBranch        *b_ph_iso855_n;   //!
   TBranch        *b_ph_iso1077_n;   //!
   TBranch        *b_ph_iso1299_n;   //!
   TBranch        *b_ph_iso151111_n;   //!
   TBranch        *b_ph_iso201616_n;   //!
   TBranch        *b_ph_trigMatch_el;   //!
   TBranch        *b_ph_elMinDR;   //!
   TBranch        *b_leadPhot_pt;   //!
   TBranch        *b_sublPhot_pt;   //!
   TBranch        *b_leadPhot_lepDR;   //!
   TBranch        *b_sublPhot_lepDR;   //!
   TBranch        *b_ph_phDR;   //!
   TBranch        *b_phPhot_lepDR;   //!
   TBranch        *b_leadPhot_lepDPhi;   //!
   TBranch        *b_sublPhot_lepDPhi;   //!
   TBranch        *b_ph_phDPhi;   //!
   TBranch        *b_phPhot_lepDPhi;   //!
   TBranch        *b_dphi_met_lep1;   //!
   TBranch        *b_dphi_met_lep2;   //!
   TBranch        *b_dphi_met_ph1;   //!
   TBranch        *b_dphi_met_ph2;   //!
   TBranch        *b_mt_lep_met;   //!
   TBranch        *b_mt_lepph1_met;   //!
   TBranch        *b_mt_lepph2_met;   //!
   TBranch        *b_mt_lepphph_met;   //!
   TBranch        *b_m_leplep;   //!
   TBranch        *b_m_mumu;   //!
   TBranch        *b_m_elel;   //!
   TBranch        *b_m_leplep_uncorr;   //!
   TBranch        *b_m_lepph1;   //!
   TBranch        *b_m_lepph2;   //!
   TBranch        *b_m_lep2ph1;   //!
   TBranch        *b_m_lep2ph2;   //!
   TBranch        *b_m_lepphlead;   //!
   TBranch        *b_m_lepphsubl;   //!
   TBranch        *b_m_lep2phlead;   //!
   TBranch        *b_m_lep2phsubl;   //!
   TBranch        *b_m_leplepph;   //!
   TBranch        *b_m_leplepphph;   //!
   TBranch        *b_m_leplepph1;   //!
   TBranch        *b_m_leplepph2;   //!
   TBranch        *b_m_lepphph;   //!
   TBranch        *b_m_leplepZ;   //!
   TBranch        *b_m_3lep;   //!
   TBranch        *b_m_4lep;   //!
   TBranch        *b_pt_leplep;   //!
   TBranch        *b_pt_lepph1;   //!
   TBranch        *b_pt_lepph2;   //!
   TBranch        *b_pt_lepphph;   //!
   TBranch        *b_pt_leplepph;   //!
   TBranch        *b_pt_secondLepton;   //!
   TBranch        *b_pt_thirdLepton;   //!
   TBranch        *b_leadPhot_leadLepDR;   //!
   TBranch        *b_leadPhot_sublLepDR;   //!
   TBranch        *b_sublPhot_leadLepDR;   //!
   TBranch        *b_sublPhot_sublLepDR;   //!
   TBranch        *b_dr_ph1_leadLep;   //!
   TBranch        *b_dr_ph1_sublLep;   //!
   TBranch        *b_dr_ph2_leadLep;   //!
   TBranch        *b_dr_ph2_sublLep;   //!
   TBranch        *b_dphi_ph1_leadLep;   //!
   TBranch        *b_dphi_ph1_sublLep;   //!
   TBranch        *b_dphi_ph2_leadLep;   //!
   TBranch        *b_dphi_ph2_sublLep;   //!
   TBranch        *b_m_ph1_ph2;   //!
   TBranch        *b_dr_ph1_ph2;   //!
   TBranch        *b_dphi_ph1_ph2;   //!
   TBranch        *b_pt_ph1_ph2;   //!
   TBranch        *b_m_leadLep_ph1_ph2;   //!
   TBranch        *b_m_leadLep_ph1;   //!
   TBranch        *b_m_leadLep_ph2;   //!
   TBranch        *b_pt_leadph12;   //!
   TBranch        *b_pt_sublph12;   //!
   TBranch        *b_eta_leadph12;   //!
   TBranch        *b_eta_sublph12;   //!
   TBranch        *b_hasPixSeed_leadph12;   //!
   TBranch        *b_hasPixSeed_sublph12;   //!
   TBranch        *b_sieie_leadph12;   //!
   TBranch        *b_sieie_sublph12;   //!
   TBranch        *b_chIsoCorr_leadph12;   //!
   TBranch        *b_chIsoCorr_sublph12;   //!
   TBranch        *b_neuIsoCorr_leadph12;   //!
   TBranch        *b_neuIsoCorr_sublph12;   //!
   TBranch        *b_phoIsoCorr_leadph12;   //!
   TBranch        *b_phoIsoCorr_sublph12;   //!
   TBranch        *b_isEB_leadph12;   //!
   TBranch        *b_isEB_sublph12;   //!
   TBranch        *b_isEE_leadph12;   //!
   TBranch        *b_isEE_sublph12;   //!
   TBranch        *b_truthMatchPh_leadph12;   //!
   TBranch        *b_truthMatchPh_sublph12;   //!
   TBranch        *b_truthMatchPhMomPID_leadph12;   //!
   TBranch        *b_truthMatchPhMomPID_sublph12;   //!
   TBranch        *b_m_nearestToZ;   //!
   TBranch        *b_m_minZdifflepph;   //!
   TBranch        *b_truelep_n;   //!
   TBranch        *b_trueph_n;   //!
   TBranch        *b_trueph_wmother_n;   //!
   TBranch        *b_truegenph_n;   //!
   TBranch        *b_truegenphpt15_n;   //!
   TBranch        *b_truegenphpt15WZMom;   //!
   TBranch        *b_truegenphpt15LepMom_n;   //!
   TBranch        *b_truegenphpt15QMom_n;   //!
   TBranch        *b_truelep_pt;   //!
   TBranch        *b_truelep_eta;   //!
   TBranch        *b_truelep_phi;   //!
   TBranch        *b_truelep_e;   //!
   TBranch        *b_truelep_isElec;   //!
   TBranch        *b_truelep_isMuon;   //!
   TBranch        *b_truelep_motherPID;   //!
   TBranch        *b_trueph_pt;   //!
   TBranch        *b_trueph_eta;   //!
   TBranch        *b_trueph_phi;   //!
   TBranch        *b_trueph_motherPID;   //!
   TBranch        *b_trueph_parentage;   //!
   TBranch        *b_trueph_nearestLepDR;   //!
   TBranch        *b_trueph_nearestQrkDR;   //!
   TBranch        *b_trueW_pt;   //!
   TBranch        *b_trueW_eta;   //!
   TBranch        *b_trueW_phi;   //!
   TBranch        *b_trueW_e;   //!
   TBranch        *b_trueleadlep_pt;   //!
   TBranch        *b_truesubllep_pt;   //!
   TBranch        *b_true_m_leplep;   //!
   TBranch        *b_trueleadlep_leadPhotDR;   //!
   TBranch        *b_trueleadlep_sublPhotDR;   //!
   TBranch        *b_truesubllep_leadPhotDR;   //!
   TBranch        *b_truesubllep_sublPhotDR;   //!
   TBranch        *b_truephph_dr;   //!
   TBranch        *b_truephph_dphi;   //!
   TBranch        *b_truephph_m;   //!
   TBranch        *b_truelepphph_m;   //!
   TBranch        *b_el_trigSF;   //!
   TBranch        *b_el_trigSFUP;   //!
   TBranch        *b_el_trigSFDN;   //!
   TBranch        *b_ph_idSF;   //!
   TBranch        *b_ph_idSFUP;   //!
   TBranch        *b_ph_idSFDN;   //!
   TBranch        *b_ph_evetoSF;   //!
   TBranch        *b_ph_evetoSFUP;   //!
   TBranch        *b_ph_evetoSFDN;   //!
   TBranch        *b_mu_trigSF;   //!
   TBranch        *b_mu_trigSFUP;   //!
   TBranch        *b_mu_trigSFDN;   //!
   TBranch        *b_mu_isoSF;   //!
   TBranch        *b_mu_isoSFUP;   //!
   TBranch        *b_mu_isoSFDN;   //!
   TBranch        *b_mu_idSF;   //!
   TBranch        *b_mu_idSFUP;   //!
   TBranch        *b_mu_idSFDN;   //!
   TBranch        *b_xfx_first_NNPDF30_nlo_nf_5_pdfas;   //!
   TBranch        *b_xfx_second_NNPDF30_nlo_nf_5_pdfas;   //!
   TBranch        *b_xfx_first_CT10nlo;   //!
   TBranch        *b_xfx_second_CT10nlo;   //!
   TBranch        *b_xfx_first_MSTW2008nlo68cl;   //!
   TBranch        *b_xfx_second_MSTW2008nlo68cl;   //!

   RecoHistograms(TTree *tree=0);
   virtual ~RecoHistograms();
   virtual void     Init(TTree *tree);
   virtual void     Loop();


};

#endif

#ifdef reco_histograms_cxx
RecoHistograms::RecoHistograms(TTree *tree) : fChain(0) 
{
   Init(tree);
}

RecoHistograms::~RecoHistograms()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void RecoHistograms::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LHEWeight_weights = 0;
   LHEWeight_ids = 0;
   mcPID = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcParentage = 0;
   mcStatus = 0;
   nPU = 0;
   puTrue = 0;
   el_pt = 0;
   el_eta = 0;
   el_sceta = 0;
   el_phi = 0;
   el_e = 0;
   el_pt_uncorr = 0;
   el_e_uncorr = 0;
   el_pfiso30 = 0;
   el_pfiso40 = 0;
   el_charge = 0;
   el_triggerMatch = 0;
   el_hasMatchedConv = 0;
   el_passTight = 0;
   el_passMedium = 0;
   el_passLoose = 0;
   el_passVeryLoose = 0;
   el_passTightTrig = 0;
   el_passMvaTrig = 0;
   el_passMvaNonTrig = 0;
   el_passMvaTrigNoIso = 0;
   el_passMvaNonTrigNoIso = 0;
   el_passMvaTrigOnlyIso = 0;
   el_passMvaNonTrigOnlyIso = 0;
   el_truthMatch_el = 0;
   el_truthMinDR_el = 0;
   el_truthMatchPt_el = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_e = 0;
   mu_pt_uncorr = 0;
   mu_eta_uncorr = 0;
   mu_phi_uncorr = 0;
   mu_e_uncorr = 0;
   mu_corrIso = 0;
   mu_charge = 0;
   mu_triggerMatch = 0;
   mu_triggerMatchDiMu = 0;
   mu_passTight = 0;
   mu_passTightNoIso = 0;
   mu_passTightNoD0 = 0;
   mu_passTightNoIsoNoD0 = 0;
   mu_truthMatch = 0;
   mu_truthMinDR = 0;
   ph_pt = 0;
   ph_eta = 0;
   ph_sceta = 0;
   ph_phi = 0;
   ph_e = 0;
   ph_scE = 0;
   ph_pt_uncorr = 0;
   ph_HoverE = 0;
   ph_HoverE12 = 0;
   ph_sigmaIEIE = 0;
   ph_chIsoCorr = 0;
   ph_neuIsoCorr = 0;
   ph_phoIsoCorr = 0;
   ph_eleVeto = 0;
   ph_hasPixSeed = 0;
   ph_isConv = 0;
   ph_passTight = 0;
   ph_passMedium = 0;
   ph_passLoose = 0;
   ph_passLooseNoSIEIE = 0;
   ph_passHOverELoose = 0;
   ph_passHOverEMedium = 0;
   ph_passHOverETight = 0;
   ph_passSIEIELoose = 0;
   ph_passSIEIEMedium = 0;
   ph_passSIEIETight = 0;
   ph_passChIsoCorrLoose = 0;
   ph_passChIsoCorrMedium = 0;
   ph_passChIsoCorrTight = 0;
   ph_passNeuIsoCorrLoose = 0;
   ph_passNeuIsoCorrMedium = 0;
   ph_passNeuIsoCorrTight = 0;
   ph_passPhoIsoCorrLoose = 0;
   ph_passPhoIsoCorrMedium = 0;
   ph_passPhoIsoCorrTight = 0;
   ph_truthMatch_el = 0;
   ph_truthMinDR_el = 0;
   ph_truthMatchPt_el = 0;
   ph_truthMatch_ph = 0;
   ph_truthMinDR_ph = 0;
   ph_truthMatchPt_ph = 0;
   ph_truthMatchMotherPID_ph = 0;
   ph_IsEB = 0;
   ph_IsEE = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_JECUnc = 0;
   jet_e = 0;
   jet_PUIDLoose = 0;
   jet_PUIDMedium = 0;
   jet_PUIDTight = 0;
   jet_CSV = 0;
   jet_genIndex = 0;
   jet_genPt = 0;
   jet_genEta = 0;
   jet_genPhi = 0;
   jet_genE = 0;
   ph_hasMatchedEle = 0;
   ph_trigMatch_el = 0;
   ph_elMinDR = 0;
   truelep_pt = 0;
   truelep_eta = 0;
   truelep_phi = 0;
   truelep_e = 0;
   truelep_isElec = 0;
   truelep_isMuon = 0;
   truelep_motherPID = 0;
   trueph_pt = 0;
   trueph_eta = 0;
   trueph_phi = 0;
   trueph_motherPID = 0;
   trueph_parentage = 0;
   trueph_nearestLepDR = 0;
   trueph_nearestQrkDR = 0;
   trueW_pt = 0;
   trueW_eta = 0;
   trueW_phi = 0;
   trueW_e = 0;
   xfx_first_NNPDF30_nlo_nf_5_pdfas = 0;
   xfx_second_NNPDF30_nlo_nf_5_pdfas = 0;
   xfx_first_CT10nlo = 0;
   xfx_second_CT10nlo = 0;
   xfx_first_MSTW2008nlo68cl = 0;
   xfx_second_MSTW2008nlo68cl = 0;
   // Set branch addresses and branch pointers

   //PU Reweights
   vector<string> pu_reweight_names = CutValues::PU_REWEIGHT_NAMES();
   for(unsigned int pu_name_index = 0; pu_name_index < pu_reweight_names.size(); pu_name_index++){
     PU_Reweights[pu_reweight_names[pu_name_index]] = 0;
   }

   // PDF Reweights 
   vector<string> pdf_set_names = CutValues::PDF_REWEIGHT_NAMES();
   for(vector<string>::iterator it = pdf_set_names.begin(); it != pdf_set_names.end(); it++){
     PDF_Reweights[*it].first = 0;
     PDF_Reweights[*it].second = 0;
   }


   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("LHEWeight_weights", &LHEWeight_weights, &b_LHEWeight_weights);
   fChain->SetBranchAddress("LHEWeight_ids", &LHEWeight_ids, &b_LHEWeight_ids);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfType01MET", &pfType01MET, &b_pfType01MET);
   fChain->SetBranchAddress("pfType01METPhi", &pfType01METPhi, &b_pfType01METPhi);
   fChain->SetBranchAddress("pfType01METsumEt", &pfType01METsumEt, &b_pfType01METsumEt);
   fChain->SetBranchAddress("pfType01METmEtSig", &pfType01METmEtSig, &b_pfType01METmEtSig);
   fChain->SetBranchAddress("pfType01METSig", &pfType01METSig, &b_pfType01METSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("vtx_n", &vtx_n, &b_vtx_n);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_sceta", &el_sceta, &b_el_sceta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_e", &el_e, &b_el_e);
   fChain->SetBranchAddress("el_pt_uncorr", &el_pt_uncorr, &b_el_pt_uncorr);
   fChain->SetBranchAddress("el_e_uncorr", &el_e_uncorr, &b_el_e_uncorr);
   fChain->SetBranchAddress("el_pfiso30", &el_pfiso30, &b_el_pfiso30);
   fChain->SetBranchAddress("el_pfiso40", &el_pfiso40, &b_el_pfiso40);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_triggerMatch", &el_triggerMatch, &b_el_triggerMatch);
   fChain->SetBranchAddress("el_hasMatchedConv", &el_hasMatchedConv, &b_el_hasMatchedConv);
   fChain->SetBranchAddress("el_passTight", &el_passTight, &b_el_passTight);
   fChain->SetBranchAddress("el_passMedium", &el_passMedium, &b_el_passMedium);
   fChain->SetBranchAddress("el_passLoose", &el_passLoose, &b_el_passLoose);
   fChain->SetBranchAddress("el_passVeryLoose", &el_passVeryLoose, &b_el_passVeryLoose);
   fChain->SetBranchAddress("el_passTightTrig", &el_passTightTrig, &b_el_passTightTrig);
   fChain->SetBranchAddress("el_passMvaTrig", &el_passMvaTrig, &b_el_passMvaTrig);
   fChain->SetBranchAddress("el_passMvaNonTrig", &el_passMvaNonTrig, &b_el_passMvaNonTrig);
   fChain->SetBranchAddress("el_passMvaTrigNoIso", &el_passMvaTrigNoIso, &b_el_passMvaTrigNoIso);
   fChain->SetBranchAddress("el_passMvaNonTrigNoIso", &el_passMvaNonTrigNoIso, &b_el_passMvaNonTrigNoIso);
   fChain->SetBranchAddress("el_passMvaTrigOnlyIso", &el_passMvaTrigOnlyIso, &b_el_passMvaTrigOnlyIso);
   fChain->SetBranchAddress("el_passMvaNonTrigOnlyIso", &el_passMvaNonTrigOnlyIso, &b_el_passMvaNonTrigOnlyIso);
   fChain->SetBranchAddress("el_truthMatch_el", &el_truthMatch_el, &b_el_truthMatch_el);
   fChain->SetBranchAddress("el_truthMinDR_el", &el_truthMinDR_el, &b_el_truthMinDR_el);
   fChain->SetBranchAddress("el_truthMatchPt_el", &el_truthMatchPt_el, &b_el_truthMatchPt_el);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
   fChain->SetBranchAddress("mu_pt_uncorr", &mu_pt_uncorr, &b_mu_pt_uncorr);
   fChain->SetBranchAddress("mu_eta_uncorr", &mu_eta_uncorr, &b_mu_eta_uncorr);
   fChain->SetBranchAddress("mu_phi_uncorr", &mu_phi_uncorr, &b_mu_phi_uncorr);
   fChain->SetBranchAddress("mu_e_uncorr", &mu_e_uncorr, &b_mu_e_uncorr);
   fChain->SetBranchAddress("mu_corrIso", &mu_corrIso, &b_mu_corrIso);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_triggerMatch", &mu_triggerMatch, &b_mu_triggerMatch);
   fChain->SetBranchAddress("mu_triggerMatchDiMu", &mu_triggerMatchDiMu, &b_mu_triggerMatchDiMu);
   fChain->SetBranchAddress("mu_passTight", &mu_passTight, &b_mu_passTight);
   fChain->SetBranchAddress("mu_passTightNoIso", &mu_passTightNoIso, &b_mu_passTightNoIso);
   fChain->SetBranchAddress("mu_passTightNoD0", &mu_passTightNoD0, &b_mu_passTightNoD0);
   fChain->SetBranchAddress("mu_passTightNoIsoNoD0", &mu_passTightNoIsoNoD0, &b_mu_passTightNoIsoNoD0);
   fChain->SetBranchAddress("mu_truthMatch", &mu_truthMatch, &b_mu_truthMatch);
   fChain->SetBranchAddress("mu_truthMinDR", &mu_truthMinDR, &b_mu_truthMinDR);
   fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   fChain->SetBranchAddress("ph_sceta", &ph_sceta, &b_ph_sceta);
   fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   fChain->SetBranchAddress("ph_e", &ph_e, &b_ph_e);
   fChain->SetBranchAddress("ph_scE", &ph_scE, &b_ph_scE);
   fChain->SetBranchAddress("ph_pt_uncorr", &ph_pt_uncorr, &b_ph_pt_uncorr);
   fChain->SetBranchAddress("ph_HoverE", &ph_HoverE, &b_ph_HoverE);
   fChain->SetBranchAddress("ph_HoverE12", &ph_HoverE12, &b_ph_HoverE12);
   fChain->SetBranchAddress("ph_sigmaIEIE", &ph_sigmaIEIE, &b_ph_sigmaIEIE);
   fChain->SetBranchAddress("ph_chIsoCorr", &ph_chIsoCorr, &b_ph_chIsoCorr);
   fChain->SetBranchAddress("ph_neuIsoCorr", &ph_neuIsoCorr, &b_ph_neuIsoCorr);
   fChain->SetBranchAddress("ph_phoIsoCorr", &ph_phoIsoCorr, &b_ph_phoIsoCorr);
   fChain->SetBranchAddress("ph_eleVeto", &ph_eleVeto, &b_ph_eleVeto);
   fChain->SetBranchAddress("ph_hasPixSeed", &ph_hasPixSeed, &b_ph_hasPixSeed);
   fChain->SetBranchAddress("ph_isConv", &ph_isConv, &b_ph_isConv);
   fChain->SetBranchAddress("ph_passTight", &ph_passTight, &b_ph_passTight);
   fChain->SetBranchAddress("ph_passMedium", &ph_passMedium, &b_ph_passMedium);
   fChain->SetBranchAddress("ph_passLoose", &ph_passLoose, &b_ph_passLoose);
   fChain->SetBranchAddress("ph_passLooseNoSIEIE", &ph_passLooseNoSIEIE, &b_ph_passLooseNoSIEIE);
   fChain->SetBranchAddress("ph_passHOverELoose", &ph_passHOverELoose, &b_ph_passHOverELoose);
   fChain->SetBranchAddress("ph_passHOverEMedium", &ph_passHOverEMedium, &b_ph_passHOverEMedium);
   fChain->SetBranchAddress("ph_passHOverETight", &ph_passHOverETight, &b_ph_passHOverETight);
   fChain->SetBranchAddress("ph_passSIEIELoose", &ph_passSIEIELoose, &b_ph_passSIEIELoose);
   fChain->SetBranchAddress("ph_passSIEIEMedium", &ph_passSIEIEMedium, &b_ph_passSIEIEMedium);
   fChain->SetBranchAddress("ph_passSIEIETight", &ph_passSIEIETight, &b_ph_passSIEIETight);
   fChain->SetBranchAddress("ph_passChIsoCorrLoose", &ph_passChIsoCorrLoose, &b_ph_passChIsoCorrLoose);
   fChain->SetBranchAddress("ph_passChIsoCorrMedium", &ph_passChIsoCorrMedium, &b_ph_passChIsoCorrMedium);
   fChain->SetBranchAddress("ph_passChIsoCorrTight", &ph_passChIsoCorrTight, &b_ph_passChIsoCorrTight);
   fChain->SetBranchAddress("ph_passNeuIsoCorrLoose", &ph_passNeuIsoCorrLoose, &b_ph_passNeuIsoCorrLoose);
   fChain->SetBranchAddress("ph_passNeuIsoCorrMedium", &ph_passNeuIsoCorrMedium, &b_ph_passNeuIsoCorrMedium);
   fChain->SetBranchAddress("ph_passNeuIsoCorrTight", &ph_passNeuIsoCorrTight, &b_ph_passNeuIsoCorrTight);
   fChain->SetBranchAddress("ph_passPhoIsoCorrLoose", &ph_passPhoIsoCorrLoose, &b_ph_passPhoIsoCorrLoose);
   fChain->SetBranchAddress("ph_passPhoIsoCorrMedium", &ph_passPhoIsoCorrMedium, &b_ph_passPhoIsoCorrMedium);
   fChain->SetBranchAddress("ph_passPhoIsoCorrTight", &ph_passPhoIsoCorrTight, &b_ph_passPhoIsoCorrTight);
   fChain->SetBranchAddress("ph_truthMatch_el", &ph_truthMatch_el, &b_ph_truthMatch_el);
   fChain->SetBranchAddress("ph_truthMinDR_el", &ph_truthMinDR_el, &b_ph_truthMinDR_el);
   fChain->SetBranchAddress("ph_truthMatchPt_el", &ph_truthMatchPt_el, &b_ph_truthMatchPt_el);
   fChain->SetBranchAddress("ph_truthMatch_ph", &ph_truthMatch_ph, &b_ph_truthMatch_ph);
   fChain->SetBranchAddress("ph_truthMinDR_ph", &ph_truthMinDR_ph, &b_ph_truthMinDR_ph);
   fChain->SetBranchAddress("ph_truthMatchPt_ph", &ph_truthMatchPt_ph, &b_ph_truthMatchPt_ph);
   fChain->SetBranchAddress("ph_truthMatchMotherPID_ph", &ph_truthMatchMotherPID_ph, &b_ph_truthMatchMotherPID_ph);
   fChain->SetBranchAddress("ph_IsEB", &ph_IsEB, &b_ph_IsEB);
   fChain->SetBranchAddress("ph_IsEE", &ph_IsEE, &b_ph_IsEE);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_JECUnc", &jet_JECUnc, &b_jet_JECUnc);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_PUIDLoose", &jet_PUIDLoose, &b_jet_PUIDLoose);
   fChain->SetBranchAddress("jet_PUIDMedium", &jet_PUIDMedium, &b_jet_PUIDMedium);
   fChain->SetBranchAddress("jet_PUIDTight", &jet_PUIDTight, &b_jet_PUIDTight);
   fChain->SetBranchAddress("jet_CSV", &jet_CSV, &b_jet_CSV);
   fChain->SetBranchAddress("jet_genIndex", &jet_genIndex, &b_jet_genIndex);
   fChain->SetBranchAddress("jet_genPt", &jet_genPt, &b_jet_genPt);
   fChain->SetBranchAddress("jet_genEta", &jet_genEta, &b_jet_genEta);
   fChain->SetBranchAddress("jet_genPhi", &jet_genPhi, &b_jet_genPhi);
   fChain->SetBranchAddress("jet_genE", &jet_genE, &b_jet_genE);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("PUWeightDN5", &PUWeightDN5, &b_PUWeightDN5);
   fChain->SetBranchAddress("PUWeightDN10", &PUWeightDN10, &b_PUWeightDN10);
   fChain->SetBranchAddress("PUWeightUP5", &PUWeightUP5, &b_PUWeightUP5);
   fChain->SetBranchAddress("PUWeightUP6", &PUWeightUP6, &b_PUWeightUP6);
   fChain->SetBranchAddress("PUWeightUP7", &PUWeightUP7, &b_PUWeightUP7);
   fChain->SetBranchAddress("PUWeightUP8", &PUWeightUP8, &b_PUWeightUP8);
   fChain->SetBranchAddress("PUWeightUP9", &PUWeightUP9, &b_PUWeightUP9);
   fChain->SetBranchAddress("PUWeightUP10", &PUWeightUP10, &b_PUWeightUP10);
   fChain->SetBranchAddress("PUWeightUP11", &PUWeightUP11, &b_PUWeightUP11);
   fChain->SetBranchAddress("PUWeightUP12", &PUWeightUP12, &b_PUWeightUP12);
   fChain->SetBranchAddress("PUWeightUP13", &PUWeightUP13, &b_PUWeightUP13);
   fChain->SetBranchAddress("PUWeightUP14", &PUWeightUP14, &b_PUWeightUP14);
   fChain->SetBranchAddress("PUWeightUP15", &PUWeightUP15, &b_PUWeightUP15);
   fChain->SetBranchAddress("PUWeightUP16", &PUWeightUP16, &b_PUWeightUP16);
   fChain->SetBranchAddress("PUWeightUP17", &PUWeightUP17, &b_PUWeightUP17);
   fChain->SetBranchAddress("passTrig_ele27WP80", &passTrig_ele27WP80, &b_passTrig_ele27WP80);
   fChain->SetBranchAddress("passTrig_mu24eta2p1", &passTrig_mu24eta2p1, &b_passTrig_mu24eta2p1);
   fChain->SetBranchAddress("passTrig_mu24", &passTrig_mu24, &b_passTrig_mu24);
   fChain->SetBranchAddress("passTrig_mu17_mu8", &passTrig_mu17_mu8, &b_passTrig_mu17_mu8);
   fChain->SetBranchAddress("passTrig_mu17_Tkmu8", &passTrig_mu17_Tkmu8, &b_passTrig_mu17_Tkmu8);
   fChain->SetBranchAddress("passTrig_ele17_ele8_9", &passTrig_ele17_ele8_9, &b_passTrig_ele17_ele8_9);
   fChain->SetBranchAddress("passTrig_ele17_ele8_22", &passTrig_ele17_ele8_22, &b_passTrig_ele17_ele8_22);
   fChain->SetBranchAddress("isBlinded", &isBlinded, &b_isBlinded);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
   fChain->SetBranchAddress("ph_hasMatchedEle", &ph_hasMatchedEle, &b_ph_hasMatchedEle);
   fChain->SetBranchAddress("mu_pt25_n", &mu_pt25_n, &b_mu_pt25_n);
   fChain->SetBranchAddress("mu_passtrig_n", &mu_passtrig_n, &b_mu_passtrig_n);
   fChain->SetBranchAddress("mu_passtrig25_n", &mu_passtrig25_n, &b_mu_passtrig25_n);
   fChain->SetBranchAddress("el_pt25_n", &el_pt25_n, &b_el_pt25_n);
   fChain->SetBranchAddress("el_passtrig_n", &el_passtrig_n, &b_el_passtrig_n);
   fChain->SetBranchAddress("el_passtrig28_n", &el_passtrig28_n, &b_el_passtrig28_n);
   fChain->SetBranchAddress("ph_mediumNoSIEIE_n", &ph_mediumNoSIEIE_n, &b_ph_mediumNoSIEIE_n);
   fChain->SetBranchAddress("ph_medium_n", &ph_medium_n, &b_ph_medium_n);
   fChain->SetBranchAddress("ph_mediumNoEleVeto_n", &ph_mediumNoEleVeto_n, &b_ph_mediumNoEleVeto_n);
   fChain->SetBranchAddress("ph_mediumNoSIEIENoEleVeto_n", &ph_mediumNoSIEIENoEleVeto_n, &b_ph_mediumNoSIEIENoEleVeto_n);
   fChain->SetBranchAddress("ph_mediumNoChIsoNoEleVeto_n", &ph_mediumNoChIsoNoEleVeto_n, &b_ph_mediumNoChIsoNoEleVeto_n);
   fChain->SetBranchAddress("ph_mediumNoNeuIsoNoEleVeto_n", &ph_mediumNoNeuIsoNoEleVeto_n, &b_ph_mediumNoNeuIsoNoEleVeto_n);
   fChain->SetBranchAddress("ph_mediumNoPhoIsoNoEleVeto_n", &ph_mediumNoPhoIsoNoEleVeto_n, &b_ph_mediumNoPhoIsoNoEleVeto_n);
   fChain->SetBranchAddress("ph_mediumNoIso_n", &ph_mediumNoIso_n, &b_ph_mediumNoIso_n);
   fChain->SetBranchAddress("ph_mediumNoChIso_n", &ph_mediumNoChIso_n, &b_ph_mediumNoChIso_n);
   fChain->SetBranchAddress("ph_mediumNoNeuIso_n", &ph_mediumNoNeuIso_n, &b_ph_mediumNoNeuIso_n);
   fChain->SetBranchAddress("ph_mediumNoPhoIso_n", &ph_mediumNoPhoIso_n, &b_ph_mediumNoPhoIso_n);
   fChain->SetBranchAddress("ph_mediumNoChIsoNoNeuIso_n", &ph_mediumNoChIsoNoNeuIso_n, &b_ph_mediumNoChIsoNoNeuIso_n);
   fChain->SetBranchAddress("ph_mediumNoChIsoNoPhoIso_n", &ph_mediumNoChIsoNoPhoIso_n, &b_ph_mediumNoChIsoNoPhoIso_n);
   fChain->SetBranchAddress("ph_mediumNoNeuIsoNoPhoIso_n", &ph_mediumNoNeuIsoNoPhoIso_n, &b_ph_mediumNoNeuIsoNoPhoIso_n);
   fChain->SetBranchAddress("ph_iso533_n", &ph_iso533_n, &b_ph_iso533_n);
   fChain->SetBranchAddress("ph_iso855_n", &ph_iso855_n, &b_ph_iso855_n);
   fChain->SetBranchAddress("ph_iso1077_n", &ph_iso1077_n, &b_ph_iso1077_n);
   fChain->SetBranchAddress("ph_iso1299_n", &ph_iso1299_n, &b_ph_iso1299_n);
   fChain->SetBranchAddress("ph_iso151111_n", &ph_iso151111_n, &b_ph_iso151111_n);
   fChain->SetBranchAddress("ph_iso201616_n", &ph_iso201616_n, &b_ph_iso201616_n);
   fChain->SetBranchAddress("ph_trigMatch_el", &ph_trigMatch_el, &b_ph_trigMatch_el);
   fChain->SetBranchAddress("ph_elMinDR", &ph_elMinDR, &b_ph_elMinDR);
   fChain->SetBranchAddress("leadPhot_pt", &leadPhot_pt, &b_leadPhot_pt);
   fChain->SetBranchAddress("sublPhot_pt", &sublPhot_pt, &b_sublPhot_pt);
   fChain->SetBranchAddress("leadPhot_lepDR", &leadPhot_lepDR, &b_leadPhot_lepDR);
   fChain->SetBranchAddress("sublPhot_lepDR", &sublPhot_lepDR, &b_sublPhot_lepDR);
   fChain->SetBranchAddress("ph_phDR", &ph_phDR, &b_ph_phDR);
   fChain->SetBranchAddress("phPhot_lepDR", &phPhot_lepDR, &b_phPhot_lepDR);
   fChain->SetBranchAddress("leadPhot_lepDPhi", &leadPhot_lepDPhi, &b_leadPhot_lepDPhi);
   fChain->SetBranchAddress("sublPhot_lepDPhi", &sublPhot_lepDPhi, &b_sublPhot_lepDPhi);
   fChain->SetBranchAddress("ph_phDPhi", &ph_phDPhi, &b_ph_phDPhi);
   fChain->SetBranchAddress("phPhot_lepDPhi", &phPhot_lepDPhi, &b_phPhot_lepDPhi);
   fChain->SetBranchAddress("dphi_met_lep1", &dphi_met_lep1, &b_dphi_met_lep1);
   fChain->SetBranchAddress("dphi_met_lep2", &dphi_met_lep2, &b_dphi_met_lep2);
   fChain->SetBranchAddress("dphi_met_ph1", &dphi_met_ph1, &b_dphi_met_ph1);
   fChain->SetBranchAddress("dphi_met_ph2", &dphi_met_ph2, &b_dphi_met_ph2);
   fChain->SetBranchAddress("mt_lep_met", &mt_lep_met, &b_mt_lep_met);
   fChain->SetBranchAddress("mt_lepph1_met", &mt_lepph1_met, &b_mt_lepph1_met);
   fChain->SetBranchAddress("mt_lepph2_met", &mt_lepph2_met, &b_mt_lepph2_met);
   fChain->SetBranchAddress("mt_lepphph_met", &mt_lepphph_met, &b_mt_lepphph_met);
   fChain->SetBranchAddress("m_leplep", &m_leplep, &b_m_leplep);
   fChain->SetBranchAddress("m_mumu", &m_mumu, &b_m_mumu);
   fChain->SetBranchAddress("m_elel", &m_elel, &b_m_elel);
   fChain->SetBranchAddress("m_leplep_uncorr", &m_leplep_uncorr, &b_m_leplep_uncorr);
   fChain->SetBranchAddress("m_lepph1", &m_lepph1, &b_m_lepph1);
   fChain->SetBranchAddress("m_lepph2", &m_lepph2, &b_m_lepph2);
   fChain->SetBranchAddress("m_lep2ph1", &m_lep2ph1, &b_m_lep2ph1);
   fChain->SetBranchAddress("m_lep2ph2", &m_lep2ph2, &b_m_lep2ph2);
   fChain->SetBranchAddress("m_lepphlead", &m_lepphlead, &b_m_lepphlead);
   fChain->SetBranchAddress("m_lepphsubl", &m_lepphsubl, &b_m_lepphsubl);
   fChain->SetBranchAddress("m_lep2phlead", &m_lep2phlead, &b_m_lep2phlead);
   fChain->SetBranchAddress("m_lep2phsubl", &m_lep2phsubl, &b_m_lep2phsubl);
   fChain->SetBranchAddress("m_leplepph", &m_leplepph, &b_m_leplepph);
   fChain->SetBranchAddress("m_leplepphph", &m_leplepphph, &b_m_leplepphph);
   fChain->SetBranchAddress("m_leplepph1", &m_leplepph1, &b_m_leplepph1);
   fChain->SetBranchAddress("m_leplepph2", &m_leplepph2, &b_m_leplepph2);
   fChain->SetBranchAddress("m_lepphph", &m_lepphph, &b_m_lepphph);
   fChain->SetBranchAddress("m_leplepZ", &m_leplepZ, &b_m_leplepZ);
   fChain->SetBranchAddress("m_3lep", &m_3lep, &b_m_3lep);
   fChain->SetBranchAddress("m_4lep", &m_4lep, &b_m_4lep);
   fChain->SetBranchAddress("pt_leplep", &pt_leplep, &b_pt_leplep);
   fChain->SetBranchAddress("pt_lepph1", &pt_lepph1, &b_pt_lepph1);
   fChain->SetBranchAddress("pt_lepph2", &pt_lepph2, &b_pt_lepph2);
   fChain->SetBranchAddress("pt_lepphph", &pt_lepphph, &b_pt_lepphph);
   fChain->SetBranchAddress("pt_leplepph", &pt_leplepph, &b_pt_leplepph);
   fChain->SetBranchAddress("pt_secondLepton", &pt_secondLepton, &b_pt_secondLepton);
   fChain->SetBranchAddress("pt_thirdLepton", &pt_thirdLepton, &b_pt_thirdLepton);
   fChain->SetBranchAddress("leadPhot_leadLepDR", &leadPhot_leadLepDR, &b_leadPhot_leadLepDR);
   fChain->SetBranchAddress("leadPhot_sublLepDR", &leadPhot_sublLepDR, &b_leadPhot_sublLepDR);
   fChain->SetBranchAddress("sublPhot_leadLepDR", &sublPhot_leadLepDR, &b_sublPhot_leadLepDR);
   fChain->SetBranchAddress("sublPhot_sublLepDR", &sublPhot_sublLepDR, &b_sublPhot_sublLepDR);
   fChain->SetBranchAddress("dr_ph1_leadLep", &dr_ph1_leadLep, &b_dr_ph1_leadLep);
   fChain->SetBranchAddress("dr_ph1_sublLep", &dr_ph1_sublLep, &b_dr_ph1_sublLep);
   fChain->SetBranchAddress("dr_ph2_leadLep", &dr_ph2_leadLep, &b_dr_ph2_leadLep);
   fChain->SetBranchAddress("dr_ph2_sublLep", &dr_ph2_sublLep, &b_dr_ph2_sublLep);
   fChain->SetBranchAddress("dphi_ph1_leadLep", &dphi_ph1_leadLep, &b_dphi_ph1_leadLep);
   fChain->SetBranchAddress("dphi_ph1_sublLep", &dphi_ph1_sublLep, &b_dphi_ph1_sublLep);
   fChain->SetBranchAddress("dphi_ph2_leadLep", &dphi_ph2_leadLep, &b_dphi_ph2_leadLep);
   fChain->SetBranchAddress("dphi_ph2_sublLep", &dphi_ph2_sublLep, &b_dphi_ph2_sublLep);
   fChain->SetBranchAddress("m_ph1_ph2", &m_ph1_ph2, &b_m_ph1_ph2);
   fChain->SetBranchAddress("dr_ph1_ph2", &dr_ph1_ph2, &b_dr_ph1_ph2);
   fChain->SetBranchAddress("dphi_ph1_ph2", &dphi_ph1_ph2, &b_dphi_ph1_ph2);
   fChain->SetBranchAddress("pt_ph1_ph2", &pt_ph1_ph2, &b_pt_ph1_ph2);
   fChain->SetBranchAddress("m_leadLep_ph1_ph2", &m_leadLep_ph1_ph2, &b_m_leadLep_ph1_ph2);
   fChain->SetBranchAddress("m_leadLep_ph1", &m_leadLep_ph1, &b_m_leadLep_ph1);
   fChain->SetBranchAddress("m_leadLep_ph2", &m_leadLep_ph2, &b_m_leadLep_ph2);
   fChain->SetBranchAddress("pt_leadph12", &pt_leadph12, &b_pt_leadph12);
   fChain->SetBranchAddress("pt_sublph12", &pt_sublph12, &b_pt_sublph12);
   fChain->SetBranchAddress("eta_leadph12", &eta_leadph12, &b_eta_leadph12);
   fChain->SetBranchAddress("eta_sublph12", &eta_sublph12, &b_eta_sublph12);
   fChain->SetBranchAddress("hasPixSeed_leadph12", &hasPixSeed_leadph12, &b_hasPixSeed_leadph12);
   fChain->SetBranchAddress("hasPixSeed_sublph12", &hasPixSeed_sublph12, &b_hasPixSeed_sublph12);
   fChain->SetBranchAddress("sieie_leadph12", &sieie_leadph12, &b_sieie_leadph12);
   fChain->SetBranchAddress("sieie_sublph12", &sieie_sublph12, &b_sieie_sublph12);
   fChain->SetBranchAddress("chIsoCorr_leadph12", &chIsoCorr_leadph12, &b_chIsoCorr_leadph12);
   fChain->SetBranchAddress("chIsoCorr_sublph12", &chIsoCorr_sublph12, &b_chIsoCorr_sublph12);
   fChain->SetBranchAddress("neuIsoCorr_leadph12", &neuIsoCorr_leadph12, &b_neuIsoCorr_leadph12);
   fChain->SetBranchAddress("neuIsoCorr_sublph12", &neuIsoCorr_sublph12, &b_neuIsoCorr_sublph12);
   fChain->SetBranchAddress("phoIsoCorr_leadph12", &phoIsoCorr_leadph12, &b_phoIsoCorr_leadph12);
   fChain->SetBranchAddress("phoIsoCorr_sublph12", &phoIsoCorr_sublph12, &b_phoIsoCorr_sublph12);
   fChain->SetBranchAddress("isEB_leadph12", &isEB_leadph12, &b_isEB_leadph12);
   fChain->SetBranchAddress("isEB_sublph12", &isEB_sublph12, &b_isEB_sublph12);
   fChain->SetBranchAddress("isEE_leadph12", &isEE_leadph12, &b_isEE_leadph12);
   fChain->SetBranchAddress("isEE_sublph12", &isEE_sublph12, &b_isEE_sublph12);
   fChain->SetBranchAddress("truthMatchPh_leadph12", &truthMatchPh_leadph12, &b_truthMatchPh_leadph12);
   fChain->SetBranchAddress("truthMatchPh_sublph12", &truthMatchPh_sublph12, &b_truthMatchPh_sublph12);
   fChain->SetBranchAddress("truthMatchPhMomPID_leadph12", &truthMatchPhMomPID_leadph12, &b_truthMatchPhMomPID_leadph12);
   fChain->SetBranchAddress("truthMatchPhMomPID_sublph12", &truthMatchPhMomPID_sublph12, &b_truthMatchPhMomPID_sublph12);
   fChain->SetBranchAddress("m_nearestToZ", &m_nearestToZ, &b_m_nearestToZ);
   fChain->SetBranchAddress("m_minZdifflepph", &m_minZdifflepph, &b_m_minZdifflepph);
   fChain->SetBranchAddress("truelep_n", &truelep_n, &b_truelep_n);
   fChain->SetBranchAddress("trueph_n", &trueph_n, &b_trueph_n);
   fChain->SetBranchAddress("trueph_wmother_n", &trueph_wmother_n, &b_trueph_wmother_n);
   fChain->SetBranchAddress("truegenph_n", &truegenph_n, &b_truegenph_n);
   fChain->SetBranchAddress("truegenphpt15_n", &truegenphpt15_n, &b_truegenphpt15_n);
   fChain->SetBranchAddress("truegenphpt15WZMom", &truegenphpt15WZMom, &b_truegenphpt15WZMom);
   fChain->SetBranchAddress("truegenphpt15LepMom_n", &truegenphpt15LepMom_n, &b_truegenphpt15LepMom_n);
   fChain->SetBranchAddress("truegenphpt15QMom_n", &truegenphpt15QMom_n, &b_truegenphpt15QMom_n);
   fChain->SetBranchAddress("truelep_pt", &truelep_pt, &b_truelep_pt);
   fChain->SetBranchAddress("truelep_eta", &truelep_eta, &b_truelep_eta);
   fChain->SetBranchAddress("truelep_phi", &truelep_phi, &b_truelep_phi);
   fChain->SetBranchAddress("truelep_e", &truelep_e, &b_truelep_e);
   fChain->SetBranchAddress("truelep_isElec", &truelep_isElec, &b_truelep_isElec);
   fChain->SetBranchAddress("truelep_isMuon", &truelep_isMuon, &b_truelep_isMuon);
   fChain->SetBranchAddress("truelep_motherPID", &truelep_motherPID, &b_truelep_motherPID);
   fChain->SetBranchAddress("trueph_pt", &trueph_pt, &b_trueph_pt);
   fChain->SetBranchAddress("trueph_eta", &trueph_eta, &b_trueph_eta);
   fChain->SetBranchAddress("trueph_phi", &trueph_phi, &b_trueph_phi);
   fChain->SetBranchAddress("trueph_motherPID", &trueph_motherPID, &b_trueph_motherPID);
   fChain->SetBranchAddress("trueph_parentage", &trueph_parentage, &b_trueph_parentage);
   fChain->SetBranchAddress("trueph_nearestLepDR", &trueph_nearestLepDR, &b_trueph_nearestLepDR);
   fChain->SetBranchAddress("trueph_nearestQrkDR", &trueph_nearestQrkDR, &b_trueph_nearestQrkDR);
   fChain->SetBranchAddress("trueW_pt", &trueW_pt, &b_trueW_pt);
   fChain->SetBranchAddress("trueW_eta", &trueW_eta, &b_trueW_eta);
   fChain->SetBranchAddress("trueW_phi", &trueW_phi, &b_trueW_phi);
   fChain->SetBranchAddress("trueW_e", &trueW_e, &b_trueW_e);
   fChain->SetBranchAddress("trueleadlep_pt", &trueleadlep_pt, &b_trueleadlep_pt);
   fChain->SetBranchAddress("truesubllep_pt", &truesubllep_pt, &b_truesubllep_pt);
   fChain->SetBranchAddress("true_m_leplep", &true_m_leplep, &b_true_m_leplep);
   fChain->SetBranchAddress("trueleadlep_leadPhotDR", &trueleadlep_leadPhotDR, &b_trueleadlep_leadPhotDR);
   fChain->SetBranchAddress("trueleadlep_sublPhotDR", &trueleadlep_sublPhotDR, &b_trueleadlep_sublPhotDR);
   fChain->SetBranchAddress("truesubllep_leadPhotDR", &truesubllep_leadPhotDR, &b_truesubllep_leadPhotDR);
   fChain->SetBranchAddress("truesubllep_sublPhotDR", &truesubllep_sublPhotDR, &b_truesubllep_sublPhotDR);
   fChain->SetBranchAddress("truephph_dr", &truephph_dr, &b_truephph_dr);
   fChain->SetBranchAddress("truephph_dphi", &truephph_dphi, &b_truephph_dphi);
   fChain->SetBranchAddress("truephph_m", &truephph_m, &b_truephph_m);
   fChain->SetBranchAddress("truelepphph_m", &truelepphph_m, &b_truelepphph_m);
   fChain->SetBranchAddress("el_trigSF", &el_trigSF, &b_el_trigSF);
   fChain->SetBranchAddress("el_trigSFUP", &el_trigSFUP, &b_el_trigSFUP);
   fChain->SetBranchAddress("el_trigSFDN", &el_trigSFDN, &b_el_trigSFDN);
   fChain->SetBranchAddress("ph_idSF", &ph_idSF, &b_ph_idSF);
   fChain->SetBranchAddress("ph_idSFUP", &ph_idSFUP, &b_ph_idSFUP);
   fChain->SetBranchAddress("ph_idSFDN", &ph_idSFDN, &b_ph_idSFDN);
   fChain->SetBranchAddress("ph_evetoSF", &ph_evetoSF, &b_ph_evetoSF);
   fChain->SetBranchAddress("ph_evetoSFUP", &ph_evetoSFUP, &b_ph_evetoSFUP);
   fChain->SetBranchAddress("ph_evetoSFDN", &ph_evetoSFDN, &b_ph_evetoSFDN);
   fChain->SetBranchAddress("mu_trigSF", &mu_trigSF, &b_mu_trigSF);
   fChain->SetBranchAddress("mu_trigSFUP", &mu_trigSFUP, &b_mu_trigSFUP);
   fChain->SetBranchAddress("mu_trigSFDN", &mu_trigSFDN, &b_mu_trigSFDN);
   fChain->SetBranchAddress("mu_isoSF", &mu_isoSF, &b_mu_isoSF);
   fChain->SetBranchAddress("mu_isoSFUP", &mu_isoSFUP, &b_mu_isoSFUP);
   fChain->SetBranchAddress("mu_isoSFDN", &mu_isoSFDN, &b_mu_isoSFDN);
   fChain->SetBranchAddress("mu_idSF", &mu_idSF, &b_mu_idSF);
   fChain->SetBranchAddress("mu_idSFUP", &mu_idSFUP, &b_mu_idSFUP);
   fChain->SetBranchAddress("mu_idSFDN", &mu_idSFDN, &b_mu_idSFDN);
   fChain->SetBranchAddress("xfx_first_NNPDF30_nlo_nf_5_pdfas", &xfx_first_NNPDF30_nlo_nf_5_pdfas, &b_xfx_first_NNPDF30_nlo_nf_5_pdfas);
   fChain->SetBranchAddress("xfx_second_NNPDF30_nlo_nf_5_pdfas", &xfx_second_NNPDF30_nlo_nf_5_pdfas, &b_xfx_second_NNPDF30_nlo_nf_5_pdfas);
   fChain->SetBranchAddress("xfx_first_CT10nlo", &xfx_first_CT10nlo, &b_xfx_first_CT10nlo);
   fChain->SetBranchAddress("xfx_second_CT10nlo", &xfx_second_CT10nlo, &b_xfx_second_CT10nlo);
   fChain->SetBranchAddress("xfx_first_MSTW2008nlo68cl", &xfx_first_MSTW2008nlo68cl, &b_xfx_first_MSTW2008nlo68cl);
   fChain->SetBranchAddress("xfx_second_MSTW2008nlo68cl", &xfx_second_MSTW2008nlo68cl, &b_xfx_second_MSTW2008nlo68cl);

   //PU Reweights                                                                                                   
   for(unsigned int pu_name_index = 0; pu_name_index < pu_reweight_names.size(); pu_name_index++){
     string pu_reweight_name = pu_reweight_names[pu_name_index];
     fChain->SetBranchAddress(pu_reweight_name.c_str(), &PU_Reweights[pu_reweight_name]);
   }

   //PDF Reweights                                                                                                  
   for(unsigned int pdf_name_index = 0; pdf_name_index < pdf_set_names.size(); pdf_name_index++){
     string pdf_set_name = pdf_set_names[pdf_name_index];
     string pdf_first_branch_name = "xfx_first_" + pdf_set_name;
     fChain->SetBranchAddress(pdf_first_branch_name.c_str(), &(PDF_Reweights[pdf_set_name].first));
     string pdf_second_branch_name = "xfx_second_" + pdf_set_name;
     fChain->SetBranchAddress(pdf_second_branch_name.c_str(), &(PDF_Reweights[pdf_set_name].second));
   }
}

#endif // #ifdef RecoHistograms_cxx
