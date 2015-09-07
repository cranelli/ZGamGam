/*
 * Example Execution
 * ./GenHistograms.exe  "/data/users/cranelli/WGamGam/Acceptances/Dressed/CommonFiducial_NLO_wMT_Dress500MeV_Skim_PUWeights_PDFReweights/job_NLO_WAA_FSR/tree.root" "Histograms/FSR_Dressed_GEN_Histograms.root"
 */
#define gen_histograms_cxx
#include "gen_histograms.h"

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

#include <iostream>
#include <string>
#include <cmath>
#include <map>

using namespace std;
string TREE_LOC = "EventTree";

int main(int argc, char * argv[]){
  //Catch if the user does not supply enough command line arguments                                                                                   
  if(argc !=3){
    cout << "useage: " << argv[0] << " <infilename> <outfilename>" << endl;
    return 0;
  }
  string in_file_name = argv[1];
  string out_file_name = argv[2];
  
  TFile * infile = new TFile(in_file_name.c_str(), "READ");
  TTree * tree = (TTree *) infile->Get(TREE_LOC.c_str());
  
  TFile * outfile = new TFile(out_file_name.c_str(), "RECREATE");
  
  // TDirectory * skim_directory = outfile->mkdir("ggNtuplizer");
  //skim_directory->cd();                                                                                                                             
  //HistogramBuilder histograms;
  
  GenHistograms gen_histograms(tree);
  gen_histograms.Loop();
  
  gen_histograms.histogram_builder_.Write();

}


void GenHistograms::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 10000;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     fChain->GetEntry(jentry);
     if(jentry%1000 == 0) cout << jentry << endl;
     //Long64_t ientry = LoadTree(jentry);

     /*
      * PileUpReweights
      */
     
     for(map<string, Float_t>::iterator it = PU_Reweights.begin(); it != PU_Reweights.end(); it++){
       //cout << it->second << endl;
     }
     

     /*                                                                         
      * NLO Weight                                                              
      */
     double nlo_weight = LHEWeight_weights->at(0);
     
     /*                                                                                                 
      * Assign Particles                                                                                
      */
     vector<MCParticleData> photons = AssignParticleByIDandStatus(CutValues::PHOTON_PDGID,
								  CutValues::FINAL_STATE_STATUS);
     vector<MCParticleData> electrons = AssignParticleByIDandStatus(CutValues::ELECTRON_PDGID,
								    CutValues::FINAL_STATE_STATUS);
     vector<MCParticleData> muons = AssignParticleByIDandStatus(CutValues::MUON_PDGID,
								CutValues::FINAL_STATE_STATUS);
     vector<MCParticleData> neutrinos = AssignParticleByIDandStatus(CutValues::NEUTRINO_PDGIDS(),
								    CutValues::FINAL_STATE_STATUS);
     //vector<MCParticleData> taus = AssignParticleByIDandStatus(CutValues::TAU_PDGID, 2);              
     //vector<MCParticleData> ws = AssignParticleByIDandStatus(24, 1);
     
     /*                                                                                                 
      * Dressing                                                                                        
      */
     
     vector<MCParticleData> dressing_photons = ObjectCuts::SelectDressPhotons(photons);
     // Dresses Leptons, and removes photons used in dressing                                           
     vector<MCParticleData> dressed_electrons = Dress(electrons, dressing_photons);
     vector<MCParticleData> dressed_muons = Dress(muons, dressing_photons);
     
     /*                                                                                                 
      * Object Cuts                                                                                     
      */

     // Input is sample of photons, with photons used for dressing removed.                             
     vector<MCParticleData> candidate_photons = ObjectCuts::SelectPhotons(dressing_photons);
     // Candidate Leptons are selected from "dressed" leptons collections                               
     vector<MCParticleData> candidate_electrons = ObjectCuts::SelectLeptons(dressed_electrons);
     vector<MCParticleData> candidate_muons = ObjectCuts::SelectLeptons(dressed_muons);
     vector<MCParticleData> candidate_neutrinos = ObjectCuts::SelectNeutrinos(neutrinos);
     
     // Event Cuts (Already Handled by the Root File Skim)
     
     /*                                                                                                 
      * Make Histograms                                                                                 
      */
     
     string decay_type = "All";
     MakeBasicHistograms(decay_type, candidate_photons);
     
     if(CutValues::DO_UNWEIGHTED){
       MakeUnweightedHistograms(decay_type, candidate_photons);
     }
   } 
}
// Get Pileup Reweights
/*
vector<WeightNamePair> GenHistograms::GetPileUpReweights(){
  vector<WeightNamePair> pileup_reweights;
  WeightNamePair pileup_pair;
  pileup_pair.name = "PUWeight5";
  pileup_pair.weight = & PUWeightUP5;
  pileup_reweights.push_back(pileup_pair);
  return pileup_reweights;
}
*/

// Basic Histograms (used for the Acceptances no reweighting for systematics).
void GenHistograms::MakeBasicHistograms(string decay_type, vector<MCParticleData> & candidate_photons){
  double nlo_weight = LHEWeight_weights->at(0);  
  double weight = nlo_weight * PUWeight;
  string prefix = decay_type;
  
  MakeHistograms(prefix, candidate_photons, weight);
}

// Makes Histograms, all events are given a weight of 1.
void GenHistograms::MakeUnweightedHistograms(string decay_type, vector<MCParticleData> & photons){
  double weight =1;
  string prefix = decay_type + "_unweighted";
  MakeHistograms(prefix, photons, weight);
}

// Make Histograms
void GenHistograms::MakeHistograms(string prefix, vector<MCParticleData> & photons, float weight){
  MCParticleData lead_photon = SelectLead(photons);
  float lead_photon_pt = lead_photon.GetFourVector().Pt();
  histogram_builder_.FillCountHistograms(prefix, weight);
  histogram_builder_.FillPtHistograms(prefix, lead_photon_pt , weight);
  histogram_builder_.FillPtCategoryHistograms(prefix, lead_photon_pt, weight);
}

/*                            
 * Return the photon in the event, with the Lead                                             
 * Transverse Momentum
 */
MCParticleData GenHistograms::SelectLead(vector<MCParticleData> & photons){
  double maxPt = 0;
  MCParticleData lead;// = NULL; 
  for(vector<MCParticleData>::iterator it = photons.begin(); it != photons.end(); ++it){
    float pt = it->GetFourVector().Pt();
    if( pt > maxPt){
      maxPt = pt;
      lead = *it;
    }
  }
  return lead;
}

/*                        
 * Dresses leptons with photons within dR cut,                                                       
 * and creates a new four vector.
 * removes photons if they are matches to a lepton.
 * leptons must be from a W or tau parent (including missing W's)
 */

vector<MCParticleData> GenHistograms::Dress(vector<MCParticleData> & leptons,
						     vector<MCParticleData> & dressing_photons){
  vector<MCParticleData> dressed_leptons;
  for(unsigned int lepton_index =0; lepton_index < leptons.size(); lepton_index++){

    MCParticleData lepton = leptons[lepton_index];
    // Must have correct parentage
    if (!ObjectCuts::PassParentCut(lepton, CutValues::LEPTON_CANDIDATE_PARENT_PDGIDS())) continue;
    if (!ObjectCuts::PassPromptCut(lepton, CutValues::NONPROMPT_BIT_MASK)) continue;
    // Add Four Vectors and Remove photon if within deltaR cut. (have to be careful with
    // incrementing indices.)
    TLorentzVector dressed_lepton_four_vector = lepton.GetFourVector();
    //unsigned int photon_index = 0;
    //while(photon_index < dressing_photons.size()){ 
    for(vector<MCParticleData>::iterator dressing_photon = dressing_photons.begin();
        dressing_photon != dressing_photons.end();){  // ++ is in body                              
      if(lepton.GetFourVector().DeltaR(dressing_photon->GetFourVector()) < CutValues::DRESSING_DR){
        dressed_lepton_four_vector += dressing_photon->GetFourVector();

        //Remove Dressing Photon from Collection                                                    
        dressing_photon = dressing_photons.erase(dressing_photon);
      } else {
        ++dressing_photon;
      }
    }
    MCParticleData dressed_lepton = lepton;
    dressed_lepton.SetFourVector(dressed_lepton_four_vector);
    dressed_leptons.push_back(dressed_lepton);
  }
  return dressed_leptons;
}

/* 
 * Return a vector of all particles, with the designated
  * ID and Status
 */

vector<MCParticleData> GenHistograms::AssignParticleByIDandStatus(int pdgID, int status){

  vector<MCParticleData> particles;
  for(int mc_index =0; mc_index < nMC; mc_index++){
    //Make Particle if it Passes Condition
    if(abs(mcPID->at(mc_index)) == pdgID && mcStatus->at(mc_index) == status){
      MCParticleData particle = MakeParticle(mc_index);
      particles.push_back(particle);
    }
  }
  return particles;
}


// Wrapper, when there is a range of pdgID values
vector<MCParticleData> GenHistograms::AssignParticleByIDandStatus(vector <int> pdgIDs, int status){
  vector<MCParticleData> total_particles;

  for(unsigned int index=0; index < pdgIDs.size(); index++){
    vector <MCParticleData> particles = AssignParticleByIDandStatus(pdgIDs[index],status);
    //Concat on vector for each pdgID
    total_particles.insert(total_particles.end(), particles.begin(), particles.end());
  }

  return total_particles;
}

MCParticleData GenHistograms::MakeParticle(int mc_index){

  MCParticleData particle;
  TLorentzVector four_vector;

  four_vector.SetPtEtaPhiE(mcPt->at(mc_index), mcEta->at(mc_index),
                           mcPhi->at(mc_index),mcE->at(mc_index));
  particle.SetFourVector(four_vector);
  particle.SetPID(mcPID->at(mc_index));
  particle.SetMomPID(mcMomPID->at(mc_index));
  particle.SetStatus(mcStatus->at(mc_index));
  particle.SetMCParentage(mcParentage->at(mc_index));

  return particle;
}
