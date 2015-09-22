#ifndef EVENT_CUTS_H
#define EVENT_CUTS_H

#include "mc_particle_data.h"
#include <vector>
#include <iostream>

using namespace std;

class EventCuts {
 public:
  
  static bool PassLeptonMultiplicity(vector<MCParticleData> & candidate_electrons, 
				     vector<MCParticleData> & candidate_muons);
  static bool PassPhotonMultiplicity(vector<MCParticleData> & candidate_photons);
  static bool PassPhotonPhotonDeltaR(vector<MCParticleData> & candidate_photons);
  static bool PassPhotonLeptonDeltaR(vector<MCParticleData> & candidate_photons, vector<MCParticleData> & candidate_leptons);  
  static bool PassLeadPt(vector<MCParticleData> & candidate_particles, float pt_cut);
  static bool PassMass(vector<MCParticleData> & candidate_particles, float mass_cut);
  
  //static bool PassMt(MCParticleData lepton, vector<MCParticleData> & neutrinos);
  

  //static void Test();
  
};


#endif
