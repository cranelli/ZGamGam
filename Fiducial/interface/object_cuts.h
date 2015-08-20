#ifndef OBJECT_CUTS_H
#define OBJECT_CUTS_H

#include "mc_particle_data.h"
#include <vector>
//#include <list>

using namespace std;

class ObjectCuts {
 public:

  static vector<MCParticleData> SelectPhotons(vector<MCParticleData> & photons);
  static vector<MCParticleData> SelectLeptons(vector<MCParticleData> & leptons);
  static vector<MCParticleData> SelectNeutrinos(vector<MCParticleData> & neutrinos);
  
  /*
  static void SelectOnPt(vector<MCParticleData> & particles, float min_pt);
  //static void SelectOnPtTest(vector<MCParticleData> & particles, float min_pt);
  static void SelectOnEta(vector<MCParticleData> & particles, float max_eta);
  static void SelectOnParent(vector<MCParticleData> & particles, int parent_id);
  static void SelectOnPromptParentage(vector<MCParticleData> & particles, int nonprompt_bit_mask);
  */

};


#endif
