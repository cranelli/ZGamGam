//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  3 10:28:02 2015 by ROOT version 5.34/30
// from TTree EventTree/EventTree
// found on file: /data/users/cranelli/WGamGam/Acceptances/Dressed/CommonFiducial_NLO_wMT_Dress500MeV_Skim_PUWeights_PDFReweights/job_NLO_WAA_ISR/tree.root
//////////////////////////////////////////////////////////

#ifndef gen_histograms_h
#define gen_histograms_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"
#include <TLorentzVector.h>

#include "mc_particle_data.h"
#include "histogram_builder.h"
#include "cut_values.h"

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>

#include <map>


// Fixed size dimensions of array or collections stored in the TTree if any.

class GenHistograms {

 private :
  string SelectDecayType(vector<MCParticleData> & candidate_electrons, vector<MCParticleData> & candidate_muons);

 public :
  // My Stuff
   HistogramBuilder histogram_builder_;
   vector<MCParticleData> AssignParticleByIDandStatus(int pdgID, int status);
   vector<MCParticleData> AssignParticleByIDandStatus(vector<int> pdgIDs, int status);
   vector<MCParticleData> Dress(vector<MCParticleData> & leptons,
                                vector<MCParticleData> & dressing_photons);

   void MakeHistograms(string prefix, vector<MCParticleData> & photons, float weight);

   void MakeBasicHistograms(string decay_type, vector<MCParticleData> & candidate_photons);
   void MakeUnweightedHistograms(string decay_type, vector<MCParticleData> & candidate_photons);
   void MakePileUpReweightHistograms(string decay_type, vector<MCParticleData> & photons);
   void MakeNLOReweightHistograms(string decay_type, vector<MCParticleData> & photons);
   void MakeCentralPDFReweightHistograms(string decay_type, vector<MCParticleData> & photons);
   void MakeEigenvectorPDFReweightHistograms(string decay_type, vector<MCParticleData> & photons);
   double CalcPDFReweight(pair<vector<double> *, vector<double> * > orig_pdf_pair, 
			  pair<vector<double> *, vector<double> * > new_pdf_pair, int new_eigenvector_index);

   MCParticleData SelectLead(vector<MCParticleData> & photons);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *LHEWeight_weights;
   vector<string>  *LHEWeight_ids;
   Int_t           run;
   Long64_t        event;
   Float_t         pdf[7];
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx_x;
   vector<float>   *mcVtx_y;
   vector<float>   *mcVtx_z;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<int>     *mcIndex;
   vector<int>     *mcDecayType;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   vector<float>   *muPFIsoR04_PU;
   vector<float>   *muPFIsoR03_PU;
   Float_t         PUWeightUP5;
   Float_t         PUWeightUP10;
   Float_t         PUWeightDN5;
   Float_t         PUWeightDN10;
   Float_t         PUWeight;
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
   TBranch        *b_pdf;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_muPFIsoR04_PU;   //!
   TBranch        *b_muPFIsoR03_PU;   //!
   TBranch        *b_PUWeightUP5;   //!
   TBranch        *b_PUWeightUP10;   //!
   TBranch        *b_PUWeightDN5;   //!
   TBranch        *b_PUWeightDN10;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_xfx_first_NNPDF30_nlo_nf_5_pdfas;   //!
   TBranch        *b_xfx_second_NNPDF30_nlo_nf_5_pdfas;   //!
   TBranch        *b_xfx_first_CT10nlo;   //!
   TBranch        *b_xfx_second_CT10nlo;   //!
   TBranch        *b_xfx_first_MSTW2008nlo68cl;   //!
   TBranch        *b_xfx_second_MSTW2008nlo68cl;   //!

   //TBranch        *b_PUWeights;

   GenHistograms(TTree *tree);
   virtual ~GenHistograms();
   virtual void     Init(TTree *tree);
   virtual void     Loop();

   

   /*
   void MakeCheckHistograms(HistogramBuilder & histograms, string prefix,
                            vector<MCParticleData> photons,
                            vector<MCParticleData> electrons,
                            vector<MCParticleData> muons,
                            vector<MCParticleData> neutrinos);
   */

   MCParticleData MakeParticle(int mc_index);

};

#endif

#ifdef gen_histograms_cxx
GenHistograms::GenHistograms(TTree *tree) : fChain(0) 
{
   Init(tree);
}

GenHistograms::~GenHistograms()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void GenHistograms::Init(TTree *tree)
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
   mcVtx_x = 0;
   mcVtx_y = 0;
   mcVtx_z = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcDecayType = 0;
   mcParentage = 0;
   mcStatus = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   muPFIsoR04_PU = 0;
   muPFIsoR03_PU = 0;
   xfx_first_NNPDF30_nlo_nf_5_pdfas = 0;
   xfx_second_NNPDF30_nlo_nf_5_pdfas = 0;
   xfx_first_CT10nlo = 0;
   xfx_second_CT10nlo = 0;
   xfx_first_MSTW2008nlo68cl = 0;
   xfx_second_MSTW2008nlo68cl = 0;

   //PU Reweights
   vector<string> pu_reweight_names = CutValues::PU_REWEIGHT_NAMES();
   for(unsigned int pu_name_index = 0; pu_name_index < pu_reweight_names.size(); pu_name_index++){
     PU_Reweights[pu_reweight_names[pu_name_index]] = 0;
   }

   // PDF Reweights
   vector<string> pdf_set_names = CutValues::PDF_REWEIGHT_NAMES();
   for(vector<string>::iterator it = pdf_set_names.begin(); it != pdf_set_names.end(); it++){
     //for(unsigned int pdf_name_index = 0; pdf_name_index < pdf_set_names.size(); pdf_name_index++){
     PDF_Reweights[*it].first = 0;
     PDF_Reweights[*it].second = 0;
   }



   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("LHEWeight_weights", &LHEWeight_weights, &b_LHEWeight_weights);
   fChain->SetBranchAddress("LHEWeight_ids", &LHEWeight_ids, &b_LHEWeight_ids);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx_x", &mcVtx_x, &b_mcVtx_x);
   fChain->SetBranchAddress("mcVtx_y", &mcVtx_y, &b_mcVtx_y);
   fChain->SetBranchAddress("mcVtx_z", &mcVtx_z, &b_mcVtx_z);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcDecayType", &mcDecayType, &b_mcDecayType);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("muPFIsoR04_PU", &muPFIsoR04_PU, &b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR03_PU", &muPFIsoR03_PU, &b_muPFIsoR03_PU);
   fChain->SetBranchAddress("PUWeightUP5", &PUWeightUP5, &b_PUWeightUP5);
   fChain->SetBranchAddress("PUWeightUP10", &PUWeightUP10, &b_PUWeightUP10);
   fChain->SetBranchAddress("PUWeightDN5", &PUWeightDN5, &b_PUWeightDN5);
   fChain->SetBranchAddress("PUWeightDN10", &PUWeightDN10, &b_PUWeightDN10);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
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

   //Notify();
}

#endif // #ifdef gen_histograms_cxx
