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
      * Select Channel
      */    
     string channel_type = SelectChannel();
     
     /*
      * Cut on Photons Location
      * Do Not Select Events where both photons are in the EndCap (or Other)
      */
 
     if(SelectPhotonsLocation() == "EEEE" || SelectPhotonsLocation()== "Other") continue;
     // No Photon pT cut for Zgamgam
     //if(pt_leadph12 < CutValues::MIN_LEAD_PHOTON_PT) continue;

     /*
      * Make Histograms
      */ 

     MakeBasicHistograms(channel_type);
     
     if(CutValues::DO_UNWEIGHTED){
       MakeUnweightedHistograms(channel_type);
     }
     if(CutValues::DO_PILEUP_REWEIGHT){
       MakePileUpReweightHistograms(channel_type);
     }
     if(CutValues::DO_NLO_REWEIGHT){
       MakeNLOReweightHistograms(channel_type);
     }
     if(CutValues::DO_CENTRAL_PDF_REWEIGHT){
       MakeCentralPDFReweightHistograms(channel_type);
     }
     if(CutValues::DO_EIGENVECTOR_PDF_REWEIGHT){
       MakeEigenvectorPDFReweightHistograms(channel_type);
     }
     if(CutValues::DO_SCALEFACTOR_REWEIGHT){
       MakeScaleFactorReweightHistograms(channel_type);
     }
   }
}


/*
 * Make Standard Set of Histograms.
 * The lead and sub-lead reconstructed photons are stored in the 
 * tree branches pt_leadph12 and pt_sublph12 - already declared
 * in the header file.
 */
void RecoHistograms::MakeHistograms(string prefix, double weight){
  histogram_builder_.FillCountHistograms(prefix, weight);
  histogram_builder_.FillPtHistograms(prefix, pt_leadph12, weight); 
  histogram_builder_.FillPtCategoryHistograms(prefix, pt_leadph12, weight);
}


/*
 * Basic Histograms (Used for the Acceptances no reweighting for systematics)
 */
void RecoHistograms::MakeBasicHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double scale_factor = CalcScaleFactor(channel_type);
  double weight = nlo_weight * PUWeight * scale_factor;
  //cout << weight << endl;
  string prefix = channel_type + "_ScaleFactors";
  MakeHistograms(prefix, weight);

}

/*
 * Unweighted Histograms - All events are given a weight of 1
 */
void RecoHistograms::MakeUnweightedHistograms(string channel_type){
  double weight = 1;
  string prefix = channel_type +"_unweighted";
  MakeHistograms(prefix, weight);
}

/*                                                                                                              
 * PileUpReweights default PUWeight is replaced with one of the reweighting variations.                         
 */
 void RecoHistograms::MakePileUpReweightHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double scale_factor = CalcScaleFactor(channel_type);

  for(map<string, Float_t>::iterator it = PU_Reweights.begin(); it != PU_Reweights.end(); it++){
    string pu_name = it->first;
    double pu_reweight = it->second;
    
    double weight = nlo_weight * pu_reweight * scale_factor;
    string prefix = channel_type + "_"+pu_name;
    MakeHistograms(prefix, weight);
  }
}

/*                                                                                                              
 * NLO Reweight Histograms (Factorization and Renormalization)                                                  
 */

void RecoHistograms::MakeNLOReweightHistograms(string channel_type){
  double scale_factor = CalcScaleFactor(channel_type);
  vector< pair <string, int> > nlo_reweight_names_indices = CutValues::NLO_REWEIGHT_NAMES_INDICES();
  for( vector< pair <string, int> >::iterator it = nlo_reweight_names_indices.begin(); it != nlo_reweight_names_indices.end(); it++){
    string name = it->first;
    int index = it->second;
    
    double nlo_reweight = LHEWeight_weights->at(index);
    double weight = nlo_reweight * PUWeight * scale_factor;
    string prefix = channel_type + "_" + name;
    MakeHistograms(prefix, weight);
  }
}

/*
 * Central PDF Reweight.  Multiply the reweighting factor, for the                                                                  
 * ratio of the orignal pdf set's weight compared to the new pdf set's weight.                                                      
 */

void RecoHistograms::MakeCentralPDFReweightHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double scale_factor = CalcScaleFactor(channel_type);
  
  pair<vector<double> *, vector <double> * > orig_pdf_pair = PDF_Reweights[CutValues::PDF_REWEIGHT_NAMES().at(CutValues::ORIGINAL_PDF_NAME_INDEX)];
  
  for(map<string, pair<vector <double> *, vector<double> * > >::iterator it = PDF_Reweights.begin();
      it != PDF_Reweights.end(); it++){
    
    string pdf_set_name = it->first;
    pair <vector<double> * , vector<double> * > new_pdf_pair= it->second;
    
    int new_eigenvector_index = 0; //Corresponds to Central Value                                                                   
    double pdf_reweight = CalcPDFReweight(orig_pdf_pair, new_pdf_pair, new_eigenvector_index);
    double weight = nlo_weight * PUWeight * scale_factor * pdf_reweight;
    string prefix = channel_type + "_" + pdf_set_name;
    MakeHistograms(prefix, weight);
  }
}

/*
 * Eigenvector PDF Reweights. Reweights for the orignal PDF weights, to the weights corresponding to
 * the PDF set's eigenvector values.
 */
void RecoHistograms::MakeEigenvectorPDFReweightHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double scale_factor = CalcScaleFactor(channel_type);
  
  pair<vector<double> *, vector <double> * > orig_pdf_pair = PDF_Reweights[CutValues::PDF_REWEIGHT_NAMES().at(CutValues::ORIGINAL_PDF_NAME_INDEX)];
  string eigenvector_pdf_set_name = CutValues::PDF_REWEIGHT_NAMES().at(CutValues::EIGENVECTOR_PDF_NAME_INDEX);
  pair<vector<double> *, vector <double> * > eigenvector_pdf_pair = PDF_Reweights[eigenvector_pdf_set_name];


  for(unsigned eigenvector_index = 0; eigenvector_index < eigenvector_pdf_pair.first->size(); eigenvector_index++){
    double pdf_reweight = CalcPDFReweight(orig_pdf_pair, eigenvector_pdf_pair, eigenvector_index);
    double weight = nlo_weight * PUWeight * scale_factor * pdf_reweight;
    string prefix = channel_type + "_" + eigenvector_pdf_set_name + "_" + to_string(eigenvector_index);
    MakeHistograms(prefix, weight);
  }

}

/*
 * Scale Factor Reweights
 * The factors reweighted differ depending on whether it is the 
 * Electron or Muon channel.
 */

void RecoHistograms::MakeScaleFactorReweightHistograms(string channel_type){
  double nlo_weight = LHEWeight_weights->at(0);
  double scale_factor = CalcScaleFactor(channel_type);

  vector<string> channel_sf_names;
  
  if(channel_type == "ElectronChannel")  channel_sf_names = CutValues::ELECTRON_CHANNEL_SCALEFACTOR_REWEIGHT_NAMES();
  if(channel_type == "MuonChannel") channel_sf_names = CutValues::MUON_CHANNEL_SCALEFACTOR_REWEIGHT_NAMES();

  for(vector<string>::iterator it = channel_sf_names.begin(); it != channel_sf_names.end(); it++){
    string sf_name = *it;
    ReweightTriplet scale_factor_triplet = SF_Reweights[sf_name];
    
    // Vary 1 Sigma UP
    double scalefactor_reweight_up = CalcSFReweight(scale_factor_triplet.orig, scale_factor_triplet.up);
    if(scalefactor_reweight_up == 0) cout << sf_name << endl;
    string prefix_up = channel_type + "_" + sf_name + "UP";
    double weight_up = nlo_weight * PUWeight * scale_factor * scalefactor_reweight_up;
    MakeHistograms(prefix_up, weight_up);

    // Vary 1 Sigma DN
    double scalefactor_reweight_down = CalcSFReweight(scale_factor_triplet.orig, scale_factor_triplet.down);
    if(scalefactor_reweight_down == 0) cout << sf_name << endl;
    string prefix_down = channel_type + "_" + sf_name + "DN";
    double weight_down = nlo_weight * PUWeight * scale_factor * scalefactor_reweight_down;
    MakeHistograms(prefix_down, weight_down);

  }

}


/*
 * Calculates the overall Scale Factor for the event.
 * Individual Scale Factors used depends on whether it is
 * the Electron or Muon Channel.
 */
double RecoHistograms::CalcScaleFactor(string channel_type){
  double scale_factor = 1;
  vector<string> channel_sf_names;
  if(channel_type == "ElectronChannel")  channel_sf_names = CutValues::ELECTRON_CHANNEL_SCALEFACTOR_REWEIGHT_NAMES();
  if(channel_type == "MuonChannel") channel_sf_names = CutValues::MUON_CHANNEL_SCALEFACTOR_REWEIGHT_NAMES();

  for(vector<string>::iterator it = channel_sf_names.begin(); it != channel_sf_names.end(); it++){
    string sf_name = *it;
    ReweightTriplet scale_factor_triplet = SF_Reweights[sf_name];
    scale_factor *= scale_factor_triplet.orig;
  }
  /*
  if(channel_type == "ElectronChannel"){
    scale_factor = el_trigSF * ph_idSF * ph_evetoSF;
  }
  if(channel_type == "MuonChannel"){
    scale_factor = mu_trigSF*mu_isoSF*mu_idSF*ph_idSF;
  }
    */
  
  return scale_factor;
}

/*
 * Given the amount to reweight the event by, base on the ratio
 * of the new SF to the old SF.
 */
double RecoHistograms::CalcSFReweight(double orig_sf, double new_sf){
  if(orig_sf ==0 ){
    cout << "Error with Reweighting, original sf is 0" << endl;
    return 0;
  }
  return new_sf / orig_sf;
}
	



/*                                                                                                                                  
 * Calculates the factor to reweight from the original pdf set to the new pdf set.                                                  
 * The eigenvector index of the original set is always 0, the central value.                                                        
 */
double RecoHistograms::CalcPDFReweight(pair<vector<double> * , vector<double> * > orig_pdf_pair, pair<vector<double> * , vector<double> * > new_pdf_pair, int new_eigenvector_index){
  int orig_eigenvector_index = 0; // Original PDF is always at Central Value                                                        
  
  double  orig_xfx_first = orig_pdf_pair.first->at(orig_eigenvector_index);
  double  orig_xfx_second = orig_pdf_pair.second->at(orig_eigenvector_index);
  
  double  new_xfx_first = new_pdf_pair.first->at(new_eigenvector_index);
  double  new_xfx_second = new_pdf_pair.second->at(new_eigenvector_index);
  
  double reweight = (new_xfx_first * new_xfx_second) / (orig_xfx_first * orig_xfx_second);

  return reweight;
}

/*
 * Uses the trigger and number of leptons to select whether the event
 * is in the electron or muon channels.
 */
string RecoHistograms::SelectChannel(){
  // Electron Channel

  if(el_passtrig_n> 0 &&  el_n==1 && mu_n==0){
    return "ElectronChannel";
  }

  // Muon Channel
  if(mu_passtrig25_n >0 && el_n == 0 && mu_n == 1){
    return "MuonChannel";
  }
  // Otherwise
  return "Other";
}

/*
 * Categorize the Event by the Location of the Lead
 * and Subleading Photons in the Detector's Barrel
 * or EndCap.
 */
string RecoHistograms::SelectPhotonsLocation(){
  // Both in Barrel (EBEB)
  if(isEB_leadph12 and isEB_sublph12) return "EBEB";
  //   Lead in Barrel Subl in EndCap (EBEE)
  if(isEB_leadph12 and isEE_sublph12) return "EBEE";
  // Lead in EndCap Subl in Barrel (EEEB)
  if(isEE_leadph12 and isEB_sublph12) return "EEEB";
  // Both in EndCap (EEEE)
  if(isEE_leadph12 and isEE_sublph12) return "EEEE";
  // Otherwise
  return "Other";
}
