#include <vector>
#include <string>
/*
 * Using this header file, as a configuration file,
 * to hold the cut values and other information for
 * the Common Fiducial Skims.
 */

class CutValues {
 private:

 public:

  static const int PHOTON_PDGID = 22;
  static const int ELECTRON_PDGID = 11;
  static const int MUON_PDGID = 13;
  static const int TAU_PDGID = 15;
  static const int FINAL_STATE_STATUS = 1; 

  static const int NONPROMPT_BIT_MASK = 4; 
  static const float PHOTON_CANDIDATE_MIN_PT = 15;
  static const float PHOTON_CANDIDATE_MAX_ETA = 2.5;
  static const float LEPTON_CANDIDATE_MIN_PT = 25;
  static const float LEPTON_CANDIDATE_MAX_ETA = 2.5;
  
  static const int REQ_NUM_CANDIDATE_LEPTONS=1;
  static const int REQ_NUM_CANDIDATE_PHOTONS=2;
  static const float PHOTON_PHOTON_MIN_DR = 0.4;
  static const float PHOTON_LEPTON_MIN_DR = 0.4;
  static const float MIN_MT = 40;

  //Lepton Parents are  W's and Tau's + quarks and gluons

  //Wrote Functions to Store Vectors
  static vector<int> LEPTON_CANDIDATE_PARENT_PDGIDS(){
    int lep_moms[8] = {1,2,3,4,5,15, 21, 24};
    static const vector<int> LEPTON_CANDIDATE_PARENT_PDGIDS(lep_moms,
							    lep_moms + 
							    sizeof lep_moms / sizeof lep_moms[0]);
    return LEPTON_CANDIDATE_PARENT_PDGIDS;
  }

  static vector<int> NEUTRINO_CANDIDATE_PARENT_PDGIDS(){int nu_moms[8] = {1,2,3,4,5,15, 21, 24};
    static const vector<int> NEUTRINO_CANDIDATE_PARENT_PDGIDS(nu_moms,
							    nu_moms + 
							    sizeof nu_moms / sizeof nu_moms[0]);
    return NEUTRINO_CANDIDATE_PARENT_PDGIDS;
  }

  static vector<int> NEUTRINO_PDGIDS(){
    int nu_pdgIDs[3] = {12,14,16};
    static const vector<int> NEUTRINO_PDGIDS(nu_pdgIDs,
					     nu_pdgIDs + 
					     sizeof nu_pdgIDs / sizeof nu_pdgIDs[0]);
    return NEUTRINO_PDGIDS;
  }
      
  


};
