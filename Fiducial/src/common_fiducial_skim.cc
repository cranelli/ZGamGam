#define CommonFiducialSkim_cxx
#include "common_fiducial_skim.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
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

//#include "mc_particle_identification.h"

using namespace std;
string IN_FILE_NAME ="/data/users/cranelli/WGamGam/NLO_ggNtuples/TruthSkim/job_NLO_WAA_ISR_TruthSkim.root";
string ORIG_TREE_LOC = "ggNtuplizer/EventTree";

string OUT_FILE_NAME= "ISR_CommonFiducial_NLO_wMT_SKIM.root";


void MakeParticleDataHistograms(HistogramBuilder & histograms, string prefix, vector<MCParticleData> particles);

int main(){
  cout << "hello world" << endl;
  TFile * infile = new TFile(IN_FILE_NAME.c_str(), "READ");
  TTree * orig_tree = (TTree *) infile->Get(ORIG_TREE_LOC.c_str());

  TFile * outfile = new TFile(OUT_FILE_NAME.c_str(), "RECREATE");
  TDirectory * skim_directory = outfile->mkdir("ggNtuplizer");
  skim_directory->cd();
  TTree * skim_tree = orig_tree->CloneTree(0);   
  HistogramBuilder  histograms;

  CommonFiducialSkim skimmer(orig_tree);
  skimmer.Loop(skim_tree, histograms);
  skim_tree->Write();
  
  map<string, TH1 *> hists = histograms.GetHistograms();
  for(map<string,TH1 *>::iterator it = hists.begin(); it != hists.end(); ++it){
    it->second->Write();
  }
  
  //outfile->Write();
}

void CommonFiducialSkim::Loop(TTree * skim_tree, HistogramBuilder & histograms)
{
   if (fChain == 0) return;
   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = 100000;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     fChain->GetEntry(jentry);
     if(jentry%1000 == 0) cout << jentry << endl;

     //cout << "Event: " << event << endl;
     
     bool keep_event = false;

     /*
      * NLO Weight
      * Normalized to plus minus 1
      */

     int nlo_weight =1;
     if (LHEWeight_weights->at(0) < 0) nlo_weight = -1;

     /*
      * Assign Particles
      */
              
     vector<MCParticleData> photons = AssignParticleByIDandStatus(CutValues::PHOTON_PDGID, 
								  CutValues::FINAL_STATE_STATUS);     
     vector<MCParticleData> electrons = AssignParticleByIDandStatus(CutValues::ELECTRON_PDGID, 
								    CutValues::FINAL_STATE_STATUS);
     vector<MCParticleData> muons = AssignParticleByIDandStatus(CutValues::MUON_PDGID, 
								CutValues::FINAL_STATE_STATUS);
     //vector<MCParticleData> taus = AssignParticleByIDandStatus(CutValues::TAU_PDGID, 2);
     vector<MCParticleData> neutrinos = AssignParticleByIDandStatus(CutValues::NEUTRINO_PDGIDS(), 
								    CutValues::FINAL_STATE_STATUS);
     //vector<MCParticleData> ws = AssignParticleByIDandStatus(24, 1);

     MakeCheckHistograms(histograms, "Assign", photons, 
			 electrons, muons, neutrinos);


     /*
      * Object Cuts
      */
     
     vector<MCParticleData> candidate_photons = ObjectCuts::SelectPhotons(photons);   
     vector<MCParticleData> candidate_electrons = ObjectCuts::SelectLeptons(electrons);
     vector<MCParticleData> candidate_muons = ObjectCuts::SelectLeptons(muons);
     vector<MCParticleData> candidate_neutrinos = ObjectCuts::SelectNeutrinos(neutrinos);

     MakeCheckHistograms(histograms, "Cut0", candidate_photons, 
     		 candidate_electrons, candidate_muons, candidate_neutrinos);
     
     /*
      * Event Cuts
      */
     
     if(!(EventCuts::PassLeptonMultiplicity(candidate_electrons, candidate_muons))) continue;
     
     string channel;
     vector<MCParticleData> candidate_leptons;
     
     if(candidate_electrons.size() == CutValues::REQ_NUM_CANDIDATE_LEPTONS){
       channel = "Electron_Channel";
       candidate_leptons = candidate_electrons;
     }
     if(candidate_muons.size() == CutValues::REQ_NUM_CANDIDATE_LEPTONS){
       channel = "Muon_Channel"; 
       candidate_leptons = candidate_muons;
     }
     histograms.FillCutFlowHistograms(channel, 1);

     if(!(EventCuts::PassPhotonMultiplicity(candidate_photons))) continue;     
     histograms.FillCutFlowHistograms(channel, 2);

     //DeltaR Cuts
     if(!(EventCuts::PassPhotonPhotonDeltaR(candidate_photons))) continue;
     histograms.FillCutFlowHistograms(channel, 3);
     if(!(EventCuts::PassPhotonLeptonDeltaR(candidate_photons, candidate_leptons))) continue;
     histograms.FillCutFlowHistograms(channel, 4);
     //MT Cuts
     if(!(EventCuts::PassMt(candidate_leptons[0], candidate_neutrinos))) continue;
     histograms.FillCutFlowHistograms(channel, 5);
 
     keep_event = true;
     
     if(keep_event) skim_tree->Fill();
     
   }
}


/* 
 * Return a vector of all particles, with the designated
 * ID and Status
 */

vector<MCParticleData> CommonFiducialSkim::AssignParticleByIDandStatus(int pdgID, int status){
  
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
vector<MCParticleData> CommonFiducialSkim::AssignParticleByIDandStatus(vector <int> pdgIDs, int status){
  vector<MCParticleData> total_particles;

  for(unsigned int index=0; index < pdgIDs.size(); index++){
    vector <MCParticleData> particles = AssignParticleByIDandStatus(pdgIDs[index],status);
    //Concat on vector for each pdgID
    total_particles.insert(total_particles.end(), particles.begin(), particles.end());
  }

  return total_particles;
}

MCParticleData CommonFiducialSkim::MakeParticle(int mc_index){

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

void CommonFiducialSkim::MakeCheckHistograms(HistogramBuilder & histograms, string prefix, vector<MCParticleData> photons, 
					    vector<MCParticleData> electrons, vector<MCParticleData> muons, 
					    vector<MCParticleData> neutrinos){

  MakeParticleDataHistograms(histograms, prefix+"_photons", photons);
  MakeParticleDataHistograms(histograms, prefix+"_electrons", electrons);
  MakeParticleDataHistograms(histograms, prefix+"_muons", muons);
  MakeParticleDataHistograms(histograms, prefix+"_neutrinos", neutrinos);
}

void MakeParticleDataHistograms(HistogramBuilder & histograms, string prefix, vector<MCParticleData> particles){
  for( vector<MCParticleData>::iterator particle = particles.begin(); particle != particles.end(); ++particle){
    histograms.fillPtHistograms(prefix, particle->GetFourVector().Pt());
  }
}



