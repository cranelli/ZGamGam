#ifndef cut_values_h
#define cut_values_h

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

  /*
   * Dressing
   */
  static const float DRESSING_DR = 0.1;

  /*
   * PDGIDs and Stauts
   */
  static const int PHOTON_PDGID = 22;
  static const int ELECTRON_PDGID = 11;
  static const int MUON_PDGID = 13;
  static const int TAU_PDGID = 15;
  static const int FINAL_STATE_STATUS = 1; 

  /*
   * Object Cuts
   */

  static const int NONPROMPT_BIT_MASK = 4; 
  static const float PHOTON_CANDIDATE_MIN_PT = 15;
  static const float PHOTON_CANDIDATE_MAX_ETA = 2.5;
  static const float LEPTON_CANDIDATE_MIN_PT = 25;
  static const float LEPTON_CANDIDATE_MAX_ETA = 2.5;

  /*
   * Event Cuts
   */
  
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

  static vector<int> NEUTRINO_CANDIDATE_PARENT_PDGIDS(){int nu_moms[8] = {1,2,3,4,5,15,21,24};
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



  /*
   * GEN and RECO HISTOGRAMS
   */

  // Unweighted
  const static bool DO_UNWEIGHTED = true;

  // PU Reweight Names
  const static bool DO_PILEUP_REWEIGHT = true;

  static vector<string> PU_REWEIGHT_NAMES(){
    //Declared Here
    string pu_reweight_names[] = {"PUWeightUP5", "PUWeightDN5"};
    static const vector<string> PU_REWEIGHT_NAMES(pu_reweight_names, 
						  pu_reweight_names + sizeof pu_reweight_names / sizeof pu_reweight_names[0]);
    return PU_REWEIGHT_NAMES;
  }

  //NLO Reweight Indices And Names (Factorization and Renormalization LHEWeight_weights[index])

  const static bool DO_NLO_REWEIGHT = true;
  
  static vector< pair<string, int> > NLO_REWEIGHT_NAMES_INDICES(){
    //Declared Here
    pair<string, int> nlo_reweight_names_indices[] = {make_pair("Factorization_Double", 1), make_pair("Factorization_Half", 2),
						      make_pair("Renormalization_Double", 3), make_pair("Renormalization_Half", 6)};
   
    static const vector< pair <string, int> > NLO_REWEIGHT_NAMES_INDICES(nlo_reweight_names_indices, 
								      nlo_reweight_names_indices + sizeof nlo_reweight_names_indices / sizeof nlo_reweight_names_indices[0]);
    return NLO_REWEIGHT_NAMES_INDICES;
  }


  //PDF Reweight (Central and Eigenvector)
  
  const static bool DO_CENTRAL_PDF_REWEIGHT = true;
  const static bool DO_EIGENVECTOR_PDF_REWEIGHT = true;
 
  static vector<string> PDF_REWEIGHT_NAMES(){
    //Declared Here
    string pdf_reweight_names[] = {"NNPDF30_nlo_nf_5_pdfas", "CT10nlo", "MSTW2008nlo68cl"};
    static const vector<string> PDF_REWEIGHT_NAMES(pdf_reweight_names, 
						  pdf_reweight_names + sizeof pdf_reweight_names / sizeof pdf_reweight_names[0]);
    return PDF_REWEIGHT_NAMES;
  }

  const static int ORIGINAL_PDF_NAME_INDEX = 0 ; // ie  NNPDF30_nlo_nf_5_pdfas
  
  const static int EIGENVECTOR_PDF_NAME_INDEX = 0; // ie  NNPDF30_nlo_nf_5_pdfas


  /*
   * Bin Pt Edges
   */

  //static const int NUM_PT_BINS = 4;
  //static const float PT_BIN_LOW_EDGE[NUM_PT_BINS] = {15, 25, 40, 70}; 
 

};

#endif
