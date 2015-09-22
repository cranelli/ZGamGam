//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 11 22:58:12 2015 by ROOT version 5.32/00
// from TTree EventTree/Event data
// found on file: /data/users/cranelli/WGamGam/NLO_ggNtuples/Truth/job_NLO_WAA_ISR_TruthSkim.root
//////////////////////////////////////////////////////////

#ifndef CommonFiducialDressLeptonSkim_h
#define CommonFiducialDressLeptonSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1F.h"
#include <TLorentzVector.h>

#include "mc_particle_data.h"
#include "histogram_builder.h"

// Header file for the classes stored in the TTree if any.
#include <vector>
//#include <list>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class CommonFiducialDressLeptonSkim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //int   count = 0;

   // Declaration of leaf types
   vector<float>   *LHEWeight_weights;
   vector<string>  *LHEWeight_ids;
   Int_t           run;
   Long64_t        event;
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

   // List of branches
   TBranch        *b_LHEWeight_weights;   //!
   TBranch        *b_LHEWeight_ids;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
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

   CommonFiducialDressLeptonSkim(TTree *tree=0);
   virtual ~CommonFiducialDressLeptonSkim();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TTree *skim_tree, HistogramBuilder & histograms);
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);

   //My Stuff
   //vector<MCParticleData> AssignParticleByIDandStatus(vector<int> & pdgIDs, int status);
   vector<MCParticleData> AssignParticleByIDandStatus(int pdgID, int status);
   vector<MCParticleData> AssignParticleByIDandStatus(vector<int> pdgIDs, int status);
   vector<MCParticleData> Dress(vector<MCParticleData> & leptons,
				vector<MCParticleData> & dressing_photons);
   void MakeCheckHistograms(HistogramBuilder & histograms, string prefix,  
			    vector<MCParticleData> photons, 
			    vector<MCParticleData> electrons,  
			    vector<MCParticleData> muons);
			    //vector<MCParticleData> neutrinos);
   MCParticleData MakeParticle(int mc_index); 

};

#endif

#ifdef CommonFiducialDressLeptonSkim_cxx
CommonFiducialDressLeptonSkim::CommonFiducialDressLeptonSkim(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/users/cranelli/WGamGam/NLO_ggNtuples/Truth/job_NLO_WAA_ISR_TruthSkim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/users/cranelli/WGamGam/NLO_ggNtuples/Truth/job_NLO_WAA_ISR_TruthSkim.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data/users/cranelli/WGamGam/NLO_ggNtuples/Truth/job_NLO_WAA_ISR_TruthSkim.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

CommonFiducialDressLeptonSkim::~CommonFiducialDressLeptonSkim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

/*
Int_t CommonFiducialDressLeptonSkim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
*/

 /*
Long64_t CommonFiducialDressLeptonSkim::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      //Notify();
   }
   return centry;
}
 */

void CommonFiducialDressLeptonSkim::Init(TTree *tree)
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("LHEWeight_weights", &LHEWeight_weights, &b_LHEWeight_weights);
   fChain->SetBranchAddress("LHEWeight_ids", &LHEWeight_ids, &b_LHEWeight_ids);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
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
}

#endif // #ifdef CommonFiducialDressLeptonSkim_cxx
