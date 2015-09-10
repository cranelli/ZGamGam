/*
 * Example Execution
 *./RecoHistograms.exe "/data/users/cranelli/WGamGam/ReFilterFinalNtuple/NLO_LepGammaGammaFinalElandMuUnblindAll_2015_08_01_ScaleFactors_PDFReweights/job_NLO_WAA_FSR/tree.root" "Histograms/FSR_RECO_Histograms.root"
 */

#define reco_histograms_cxx
#include "reco_histograms.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"

#include "cut_values.h"
#include "mc_particle_data.h"
#include "object_cuts.h"
#include "event_cuts.h"
#include "histogram_builder.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;
string TREE_LOC = "EventTree";

int main(int argc, char * argv[]){

  // Make Sure Correct Number of Arguments are Given
  if(argc !=3){
    cout << "useage: " << argv[0] << " <infilename> <outfilename>" << endl;
    return 0;
  }

  string in_file_name = argv[1];
  string out_file_name = argv[2];

  TFile * infile = new TFile(in_file_name.c_str(), "READ");
  TTree * tree = (TTree *) infile->Get(TREE_LOC.c_str());

  TFile * outfile = new TFile(out_file_name.c_str(), "RECREATE");

  RecoHistograms reco_histograms(tree);
  reco_histograms.Loop();
  
  reco_histograms.histogram_builder_.Write();
}


void RecoHistograms::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 10000; 
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     fChain->GetEntry(jentry);
     if(jentry%1000 == 0) cout << jentry << endl;
     //Long64_t ientry = LoadTree(jentry);

     /*                                           
      * NLO Weight 
      */
     double nlo_weight = LHEWeight_weights->at(0);
   
     string channel_type = "All";
     
     MakeBasicHistograms(channel_type);

   }
}

/*
 * Basic Histograms (Used for the Acceptances no reweighting for systematics)
 */
void RecoHistograms::MakeBasicHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double weight = nlo_weight * PUWeight;
  string prefix = channel_type;

  MakeHistograms(prefix, weight);

}

/*
 * Make Standard Histograms.
 * The lead and sub-lead reconstructed photons are stored in the 
 * tree branches pt_leadph12 and pt_sublph12 - already declared
 * in the header file.
 */
void RecoHistograms::MakeHistograms(string prefix, double weight){
  histogram_builder_.FillCountHistograms(prefix, weight);
  histogram_builder_.FillPtHistograms(prefix, pt_leadph12, weight); 
  histogram_builder_.FillPtCategoryHistograms(prefix, pt_leadph12, weight);
}
