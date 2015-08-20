//Header File
#include "event_cuts.h"
#include "TLorentzVector.h"

#include "mc_particle_data.h"
#include "cut_values.h"

#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

/*
 * Pass if the event has the right number of candidate leptons.
 * (Electron and Muon Channel are completely separate channels.)
 */

bool EventCuts::PassLeptonMultiplicity(vector<MCParticleData> & candidate_electrons,
				       vector<MCParticleData> & candidate_muons){
  bool pass = false;
  if( candidate_electrons.size() == CutValues::REQ_NUM_CANDIDATE_LEPTONS && candidate_muons.size() == 0) pass =true;
  if( candidate_electrons.size() == 0 && candidate_muons.size() == CutValues::REQ_NUM_CANDIDATE_LEPTONS) pass =true;
  return pass;
}

bool EventCuts::PassPhotonMultiplicity(vector<MCParticleData> & candidate_photons){
  return (candidate_photons.size() == CutValues::REQ_NUM_CANDIDATE_PHOTONS);
}

/*
 * Check that the Candidate Photons have the proper DeltaR
 * Separation from each other
 * (There is a double counting inefficiency)
 */
bool EventCuts::PassPhotonPhotonDeltaR(vector<MCParticleData> & candidate_photons){
  for( unsigned int index1 =0; index1 < candidate_photons.size(); index1++){
    TLorentzVector photon1 = candidate_photons.at(index1).GetFourVector();
    
    for( unsigned int index2 =0; index2 < candidate_photons.size(); index2++){
      if(index1 == index2) continue; // Do not compare to itself
      TLorentzVector photon2 = candidate_photons.at(index2).GetFourVector();

      if(photon1.DeltaR(photon2) < CutValues::PHOTON_PHOTON_MIN_DR) return false;
    }
  }
  // Only True if All Pairings Pass
  return true;
}

bool EventCuts::PassPhotonLeptonDeltaR(vector<MCParticleData> & candidate_photons,
				       vector<MCParticleData> & candidate_leptons){
  for( unsigned int photon_index =0; photon_index < candidate_photons.size(); photon_index++){
    TLorentzVector photon = candidate_photons.at(photon_index).GetFourVector();
    
    for( unsigned int lepton_index =0; lepton_index < candidate_leptons.size(); lepton_index++){
      TLorentzVector lepton = candidate_leptons.at(lepton_index).GetFourVector();
      if(photon.DeltaR(lepton) < CutValues::PHOTON_LEPTON_MIN_DR) return false;
    }
  }
  // Only True if All Pairings Pass
  return true;
}

bool EventCuts::PassMt(MCParticleData candidate_lepton, vector<MCParticleData> & candidate_neutrinos){
  TLorentzVector sum_nu_four_vector;
  for(unsigned int nu_index = 0; nu_index < candidate_neutrinos.size(); nu_index++){
    sum_nu_four_vector += candidate_neutrinos[nu_index].GetFourVector();
  }

  TLorentzVector lepton_four_vector = candidate_lepton.GetFourVector();
  float met = sum_nu_four_vector.Pt();
  
  // Calculate Transverse Mass
  float mt2 = 2*lepton_four_vector.Et() * met * (1 - cos(lepton_four_vector.DeltaPhi(sum_nu_four_vector)));
  float mt = sqrt(mt2);

  return(mt > CutValues::MIN_MT);
}								       
								       
							       

							       


/*
void EventCuts::Test(){
  cout << "EventCuts Working" << endl;
}
*/





