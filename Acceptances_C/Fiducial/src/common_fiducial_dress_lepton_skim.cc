/*
 * Example Execution
 * ./CommonFiducialDressLeptonSkim.exe "/data/users/cranelli/ZGamGam/NLO_ggNtuples/llaa_nlo_ggNtuple_part1.root" "CommonFiducialDressLepton_NLO_SKIM.root"
 */

#define CommonFiducialDressLeptonSkim_cxx
#include "common_fiducial_dress_lepton_skim.h"

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
//string IN_FILE_NAME ="/data/users/cranelli/WGamGam/NLO_ggNtuples/TruthSkim/job_NLO_WAA_ISR_PtG500MeV_TruthSkim.root";
string ORIG_TREE_LOC = "ggNtuplizer/EventTree";
//string OUT_FILE_NAME= "ISR_CommonFiducialDressLepton_NLO_wMT_SKIM.root";


void MakeParticleDataHistograms(HistogramBuilder & histograms, string prefix, vector<MCParticleData> particles);

int main(int argc, char * argv[]){
  //Catch if the user does not supply enough command line arguments 
                                                            
  if(argc !=3){
    cout << "useage: " << argv[0] << " <infilename> <outfilename>" << endl;
    return 0;
  }
  string in_file_name = argv[1];
  string out_file_name = argv[2];

  TFile * infile = new TFile(in_file_name.c_str(), "READ");
  TTree * orig_tree = (TTree *) infile->Get(ORIG_TREE_LOC.c_str());

  TFile * outfile = new TFile(out_file_name.c_str(), "RECREATE");
  TDirectory * skim_directory = outfile->mkdir("ggNtuplizer");
  skim_directory->cd();
  TTree * skim_tree = orig_tree->CloneTree(0);   
  HistogramBuilder  histograms;

  CommonFiducialDressLeptonSkim skimmer(orig_tree);
  skimmer.Loop(skim_tree, histograms);
  skim_tree->Write();
  
  histograms.Write();

}

void CommonFiducialDressLeptonSkim::Loop(TTree * skim_tree, HistogramBuilder & histograms)
{
  
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 10000;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     fChain->GetEntry(jentry);
     if(jentry%1000 == 0) cout << jentry << endl;

     //cout << "Event: " << event << endl;
     
     bool keep_event = false;

     /*
      * NLO Weight
      * Normalized to plus minus 1
      */

     double nlo_weight =1;
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
     //vector<MCParticleData> neutrinos = AssignParticleByIDandStatus(CutValues::NEUTRINO_PDGIDS(), 
     //							    CutValues::FINAL_STATE_STATUS);

     MakeCheckHistograms(histograms, "Assign", photons, 
			 electrons, muons);

     /*
      * Dressing
      */
     
     vector<MCParticleData> dressing_photons = ObjectCuts::SelectDressPhotons(photons);
     MakeCheckHistograms(histograms, "Dress0", dressing_photons, 
			 electrons, muons);
     
     // Dresses Leptons, and removes photons used in dressing
     vector<MCParticleData> dressed_electrons = Dress(electrons, dressing_photons); 
     vector<MCParticleData> dressed_muons = Dress(muons, dressing_photons);
     MakeCheckHistograms(histograms, "Dress", dressing_photons, 
			 dressed_electrons, dressed_muons);

     /*
      * Object Cuts
      */ 

     // Input is sample of photons, with photons used for dressing removed.
     vector<MCParticleData> candidate_photons = ObjectCuts::SelectPhotons(dressing_photons);
     // Candidate Leptons are selected from "dressed" leptons collections
     vector<MCParticleData> candidate_electrons = ObjectCuts::SelectLeptons(dressed_electrons);
     vector<MCParticleData> candidate_muons = ObjectCuts::SelectLeptons(dressed_muons);
     //vector<MCParticleData> candidate_neutrinos = ObjectCuts::SelectNeutrinos(neutrinos);
     

     MakeCheckHistograms(histograms, "Cut0", candidate_photons, 
     		 candidate_electrons, candidate_muons);
     
     /*
      * Event Cuts
      */
     histograms.FillCutFlowHistograms("Electron_Channel", 0, nlo_weight);
     histograms.FillCutFlowHistograms("Muon_Channel", 0, nlo_weight);

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
     histograms.FillCutFlowHistograms(channel, 1, nlo_weight);

     if(!(EventCuts::PassPhotonMultiplicity(candidate_photons))) continue;     
     histograms.FillCutFlowHistograms(channel, 2, nlo_weight);

     // Lead Lepton Pt Cut
     if(!(EventCuts::PassLeadPt(candidate_leptons, CutValues::MIN_LEAD_LEPTON_PT))) continue;
     histograms.FillCutFlowHistograms(channel, 3, nlo_weight);

     // Di-Lepton Mass Cut
     if(!(EventCuts::PassMass(candidate_leptons, CutValues::MIN_DILEPTON_MASS))) continue;
     histograms.FillCutFlowHistograms(channel, 4, nlo_weight);

     // Lepton Lepton Delta R Cuts
     if(!(EventCuts::PassDeltaR(candidate_leptons, candidate_leptons, CutValues::LEPTON_LEPTON_MIN_DR))) continue;
     histograms.FillCutFlowHistograms(channel, 6, nlo_weight);

     //Photon Photon Delta R Cuts
     if(!(EventCuts::PassDeltaR(candidate_photons, candidate_photons, CutValues::PHOTON_PHOTON_MIN_DR))) continue;
     histograms.FillCutFlowHistograms(channel, 6, nlo_weight);

     //Photon Lepton Delta R Cuts
     if(!(EventCuts::PassDeltaR(candidate_photons, candidate_leptons, CutValues::PHOTON_LEPTON_MIN_DR))) continue;
     histograms.FillCutFlowHistograms(channel, 7, nlo_weight);

     

     
     
     keep_event = true;
     
     if(keep_event) skim_tree->Fill();
     
   }
   //cout << count << endl;
}

/*
 * Dresses leptons with photons within dR cut,
 * and creates a new four vector.
 * removes photons if they are matches to a lepton.
 * leptons must be from a W or tau parent (including missing W's)
 */

vector<MCParticleData> CommonFiducialDressLeptonSkim::Dress(vector<MCParticleData> & leptons, 
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
      //MCParticleData dressing_photon = dressing_photons[photon_index];
      if(lepton.GetFourVector().DeltaR(dressing_photon->GetFourVector()) < CutValues::DRESSING_DR){
	dressed_lepton_four_vector += dressing_photon->GetFourVector();
	
	//Remove Dressing Photon from Collection
	dressing_photon = dressing_photons.erase(dressing_photon);
      }	else {
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

vector<MCParticleData> CommonFiducialDressLeptonSkim::AssignParticleByIDandStatus(int pdgID, int status){
  
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
vector<MCParticleData> CommonFiducialDressLeptonSkim::AssignParticleByIDandStatus(vector <int> pdgIDs, int status){
  vector<MCParticleData> total_particles;

  for(unsigned int index=0; index < pdgIDs.size(); index++){
    vector <MCParticleData> particles = AssignParticleByIDandStatus(pdgIDs[index],status);
    //Concat on vector for each pdgID
    total_particles.insert(total_particles.end(), particles.begin(), particles.end());
  }

  return total_particles;
}

MCParticleData CommonFiducialDressLeptonSkim::MakeParticle(int mc_index){
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

void CommonFiducialDressLeptonSkim::MakeCheckHistograms(HistogramBuilder & histograms, 
							string prefix, 
							vector<MCParticleData> photons, 
							vector<MCParticleData> electrons, 
							vector<MCParticleData> muons){

  MakeParticleDataHistograms(histograms, prefix+"_photons", photons);
  MakeParticleDataHistograms(histograms, prefix+"_electrons", electrons);
  MakeParticleDataHistograms(histograms, prefix+"_muons", muons);
  //MakeParticleDataHistograms(histograms, prefix+"_neutrinos", neutrinos);
}

void MakeParticleDataHistograms(HistogramBuilder & histograms, string prefix, vector<MCParticleData> particles){
  for( vector<MCParticleData>::iterator particle = particles.begin(); particle != particles.end(); ++particle){
    histograms.FillPtHistograms(prefix, particle->GetFourVector().Pt());
  }
}



