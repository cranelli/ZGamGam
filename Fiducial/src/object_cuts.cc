//Header File
#include "object_cuts.h"
#include "TLorentzVector.h"

#include "mc_particle_data.h"
#include "cut_values.h"

#include <vector>
#include <functional>
//#include <list>
#include <iostream>

using namespace std;

bool PassPtCut(MCParticleData & particle, float min_pt);
bool PassEtaCut(MCParticleData & particle, float max_eta); 
bool PassPromptCut(MCParticleData & particle, int prompt_bit_mask);
bool PassParentCut(MCParticleData & particle, vector<int> parent_pdgIDs);


/*
 * Definition of Functors
 */

/* 
struct IsPtLess : public binary_function<MCParticleData, float, bool> {
  bool operator() (MCParticleData particle, float min_pt) const {
    return (particle.GetFourVector().Pt() <  min_pt);
  };
};

//Select on Pt

void ObjectCuts::SelectOnPt(vector<MCParticleData> & particles, float min_pt){
  particles.erase(remove_if(
			    particles.begin(), particles.end(), bind2nd(IsPtLess(), min_pt)),
		  particles.end());  
}
*/


/*
 * Selection Criteria for Candidate Photons
 */

vector<MCParticleData> ObjectCuts::SelectPhotons(vector<MCParticleData> & photons){
  
  vector<MCParticleData> select_photons;

  //  for(vector<MCParticleData>::iterator it = particles.begin(); it != particles.end(); ++it){
  for(unsigned int index =0; index < photons.size(); index++){
    MCParticleData  photon = photons.at(index);
    
    if(!PassPtCut(photon, CutValues::PHOTON_CANDIDATE_MIN_PT)) continue;
    if(!PassEtaCut(photon, CutValues::PHOTON_CANDIDATE_MAX_ETA)) continue;
    if(!PassPromptCut(photon, CutValues::NONPROMPT_BIT_MASK)) continue;
    
    select_photons.push_back(photon);
  }
  return select_photons;
}

vector<MCParticleData> ObjectCuts::SelectLeptons(vector<MCParticleData> & leptons){
  
  vector<MCParticleData> select_leptons;

  //  for(vector<MCParticleData>::iterator it = particles.begin(); it != particles.end(); ++it){                
  for(unsigned int index =0; index < leptons.size(); index++){
    MCParticleData  lepton = leptons.at(index);

    if(!PassPtCut(lepton, CutValues::LEPTON_CANDIDATE_MIN_PT)) continue;
    if(!PassEtaCut(lepton, CutValues::LEPTON_CANDIDATE_MAX_ETA)) continue;
    if(!PassPromptCut(lepton, CutValues::NONPROMPT_BIT_MASK)) continue;
    if(!PassParentCut(lepton, CutValues::LEPTON_CANDIDATE_PARENT_PDGIDS())) continue;
       
    select_leptons.push_back(lepton);
  }
  return select_leptons;
}

vector<MCParticleData> ObjectCuts::SelectNeutrinos(vector<MCParticleData> & neutrinos){
  
  vector<MCParticleData> select_neutrinos;

  //  for(vector<MCParticleData>::iterator it = particles.begin(); it != particles.end(); ++it){
  for(unsigned int index =0; index < neutrinos.size(); index++){
    MCParticleData  neutrino = neutrinos.at(index);
    if(!PassPromptCut(neutrino, CutValues::NONPROMPT_BIT_MASK)) continue;
    if(!PassParentCut(neutrino, CutValues::NEUTRINO_CANDIDATE_PARENT_PDGIDS())) continue;
    
    select_neutrinos.push_back(neutrino);
  }
  return select_neutrinos;
}


bool PassPtCut(MCParticleData & particle, float min_pt){
  return (particle.GetFourVector().Pt() > min_pt);
}

bool PassEtaCut(MCParticleData & particle, float max_eta){
  return (abs(particle.GetFourVector().Eta()) < max_eta);
}

bool PassPromptCut(MCParticleData & particle, int non_prompt_bit_mask){
  return ((particle.GetMCParentage() & non_prompt_bit_mask) != non_prompt_bit_mask);
}

/*
 * If any of the pdgIDs match, return true.
 */
bool PassParentCut(MCParticleData & particle, vector<int> parent_pdgIDs){
  bool pass = false;
  for(unsigned int index=0; index < parent_pdgIDs.size(); index++){
    if(abs(particle.GetMomPID()) == parent_pdgIDs[index]) pass = true;
  }
  return pass;
   
}




/*
void ObjectCuts::SelectOnEta(vector<MCParticleData> & particles, float max_eta){
  particles.erase(remove_if(
			    particles.begin(), particles.end(), bind2nd(IsEtaGreater(), max_eta)),
		  particles.end());
  
}

void ObjectCuts::SelectOnParent(vector<MCParticleData> & particles, int parent_id){
  particles.erase(remove_if(
                            particles.begin(), particles.end(), bind2nd(IsMomIDUnequal(), parent_id)),
                  particles.end());
}

void ObjectCuts::SelectOnPromptParentage(vector<MCParticleData> & particles, int nonprompt_bit_mask){
  particles.erase(remove_if(
			    particles.begin(), particles.end(), bind2nd(IsNonPrompt(), nonprompt_bit_mask)),
		  particles.end());
  
}
*/

/*
class IsPtLessThan {
private:
  float min_pt_;
public:
  IsPtLessThan(float min_pt = 0.0) : min_pt_(min_pt) {}
  bool operator() (MCParticleData particle) const
  {
    return (particle.GetFourVector().Pt() <  min_pt_);
  }
};
*/

// Outdated versions using lists


// Iterating over and removing elements from list
/*
void ObjectCuts::SelectOnEta(list<MCParticleData> & particles, float max_eta){
  list<MCParticleData>::iterator it = particles.begin();
  while(it != particles.end()){
    float eta = it->GetFourList().Eta();
    if( eta > max_eta){
      particles.erase(it++);
    } 
    else {
      ++it;
    }
  }
}
*/


