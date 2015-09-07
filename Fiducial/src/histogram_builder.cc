#include "histogram_builder.h"
/*
 * The HistogramBuilder Class contains  
 * functions to build generic histograms of                           
 * different types.  The type of object (or a specific                      
 * selection cut) can be sepecified by the key.
 *                           
 * Created by Christopher Anelli
 * On 6.15.2014
*/

#include "TH1F.h"
#include "TH2F.h"
#include "math.h"

//#include "cut_values.h"

#include <string>
#include <iostream>
//#include "HoMuonTrigger/hoTriggerAnalyzer/interface/CommonFunctions.h"

using namespace std;


// Constructor
/*
HistogramBuilder::HistogramBuilder(){
};
*/



/*
 * CutFlow Histograms
 */

void HistogramBuilder::FillCutFlowHistograms(string prefix, int cut_step, double weight){
  string key = prefix + "_CutFlow";

  if(!histograms_.count(key)){
    histograms_[key] = new TH1F(Form("%s",key.c_str()),
			       Form("%s",key.c_str()),
			       10, 0, 10);
    histograms_[key]->GetYaxis()->SetTitle("CutFlow");
    histograms_[key]->Sumw2();
  }
  histograms_[key]->Fill(cut_step, weight);
}

/*
 *Pt Histograms
 */
void HistogramBuilder::FillPtHistograms(string prefix, float pt, double weight){
  string key = prefix+"_Pt";
  if(!histograms_.count(key)){
        
    histograms_[key] = new TH1F(Form("%s",key.c_str()),
				Form("%s",key.c_str()),
				100,0,200);
    SetAxises(histograms_[key],"Pt (GeV)", "Counts");
    histograms_[key]->Sumw2();
  } 
  histograms_[key]->Fill(pt, weight);
}

void HistogramBuilder::FillPtCategoryHistograms(string prefix, 
						float pt, double weight){
  string key = prefix+"_Category_Pt";
  if(!histograms_.count(key)){
        
    histograms_[key] = new TH1F(Form("%s",key.c_str()),
				Form("%s",key.c_str()),
				num_photon_pt_bins, 
				photon_pt_bin_low_edges);
    SetAxises(histograms_[key],"Pt (GeV)", "Counts");
    histograms_[key]->Sumw2();
  } 
  histograms_[key]->Fill(pt, weight);
}

void HistogramBuilder::FillEtaHistograms(string prefix, float eta, double weight){
  string key = prefix + "_Eta";

  if(!histograms_.count(key)){
    histograms_[key] = new TH1F(Form("%s",key.c_str()),
				Form("%s",key.c_str()),
				100, -3.0, 3.0);  
    SetAxises(histograms_[key], "#eta", "Counts");
    histograms_[key]->Sumw2();
    
  }
  histograms_[key]->Fill(eta, weight);
}

/*                                                                             
 *Eta Phi Histograms                                                         
 */
/*
void HistogramBuilder::fillEtaPhiHistograms(float eta, float phi, string key, double weight){


  if(!_h1Phi.count(key)){
    _h1Phi[key] = new TH1F(Form("%s_Phi",key.c_str()),
                                           Form("%s Phi",key.c_str()),
                                           500, -3.2, 3.2);  //HO has 72 iphis and 30 ietas                         
    SetAxises(_h1Phi[key], "#phi", "Counts");
    _h1Phi[key]->Sumw2();
  }
  _h1Phi[key]->Fill(phi, weight);

  if(!_h2EtaPhiMap.count(key)){
    _h2EtaPhiMap[key] = new TH2F(Form("%s_EtaPhi",key.c_str()),
                                                 Form("%s_EtaPhi",key.c_str()),
						 500, -1.5, 1.5, 500, -3.14, 3.14);
    SetAxises(_h2EtaPhiMap[key], "#eta","#phi");
    _h2EtaPhiMap[key]->Sumw2();
  }
  _h2EtaPhiMap[key]->Fill(eta, phi, weight);

};
*/

/*
 *Counting Histograms
 *Fills the 1 bin.
 */

void HistogramBuilder::FillCountHistograms(string prefix, double weight){
  string key = prefix + "_Count";
  if(!histograms_.count(key)){                                                   
    histograms_[key] = new TH1F(Form("%s_Count",key.c_str()),    
                                               Form("%s Count",key.c_str()),    
                                               2, 0, 2);
    histograms_[key]->GetYaxis()->SetTitle("Counts");
    histograms_[key]->Sumw2();
  }                                                                             
  histograms_[key]->Fill(1, weight);
}
                                                                             

/*      
 *Weight Histograms
 */
/*
void HistogramBuilder::fillWeightHistograms(float weight_val,string key, double weight){ 
  if(!_h1Weight.count(key)){                                                   
    _h1Weight[key] = new TH1F(Form("%s_Weight",key.c_str()),     
					      Form("%s Weight",key.c_str()),  
					      15000, 0, 150);
    SetAxises(_h1Weight[key],"Weight", "Counts");
  }                                                      
                       
  _h1Weight[key]->Fill(weight_val, weight);                                        
     
}
*/
     





/*
 *Delta Eta Delta Phi Histograms
 */
/*
void HistogramBuilder::fillDeltaEtaDeltaPhiHistograms(float eta1, float eta2, 
						      float phi1, float phi2, 
						      string key, double weight){
  float deltaEta, deltaPhi;
  //CommonFunctions commonFunctions;
  deltaEta = eta1 - eta2;
  //deltaPhi = commonFunctions.WrapCheck(Phi1, phi2);
  deltaPhi = WrapCheck(phi1, phi2);

  //Delta Eta Histograms Fill
  if(!_h1DeltaEta.count(key)){
        _h1DeltaEta[key] =  new TH1F(Form("%s_DeltaEta",
							  key.c_str()), 
						     Form("#Delta #eta %s",
							  key.c_str()),   
						     2000, -2.6, 2.6);
	_h1DeltaEta[key]->Sumw2();
	SetAxises(_h1DeltaEta[key],"#Delta #eta", "Counts");
  }
  _h1DeltaEta[key]->Fill(deltaEta, weight);
    
  //Delta Phi Histograms Fill
  if(!_h1DeltaPhi.count(key)){
    _h1DeltaPhi[key] = new TH1F(Form("%s_DeltaPhi",
						     key.c_str()), 
						Form("%s #Delta #Phi",
						     key.c_str()),    
                                                2000, -3.2, 3.2);
    _h1DeltaPhi[key]->Sumw2();
    SetAxises(_h1DeltaPhi[key],"#Delta #phi", "Counts");
  }
  _h1DeltaPhi[key]->Fill(deltaPhi, weight);
  
  //DeltaEta Delta Phi Histograms Fill
  if(!_h2DeltaEtaDeltaPhi.count(key)){
    _h2DeltaEtaDeltaPhi[key] = new TH2F(Form("%s_DeltaEtaDeltaPhi",key.c_str()),Form("%s #Delta#eta #Delta#Phi",key.c_str()), 2000, -2.6, 2.6, 2000, -3.14, 3.14);
    SetAxises(_h2DeltaEtaDeltaPhi[key],"#Delta #eta", "#Delta #phi");
    _h2DeltaEtaDeltaPhi[key]->Sumw2();
  }
  _h2DeltaEtaDeltaPhi[key]->Fill(deltaEta, deltaPhi, weight);
} 
*/

/*
 * Scatter Pt Plot
 */
/*
void HistogramBuilder::fillScatterPt(float pt1, float pt2, string key, double weight){
  if(!_h2ScatterPt.count(key)){
    _h2ScatterPt[key] = new TH2F(Form("%s_ScatterPt",key.c_str()),
				 Form("%s_ScatterPt",key.c_str()),
				 400, 0, 200, 400, 0, 200);
    SetAxises(_h2ScatterPt[key], "P_t (GeV)","P_t (GeV)");
    _h2ScatterPt[key]->Sumw2();
  }
  _h2ScatterPt[key]->Fill(pt1, pt2, weight);
}

*/

/*
 *
 * Helper Functions
 *
 */


void HistogramBuilder::SetAxises(TH1 * h1, string xTitle, string yTitle){
  // X Axis
  h1->GetXaxis()->SetTitle(xTitle.c_str());
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitleOffset(0.9);
  
  // Y Axis
  h1->GetYaxis()->SetTitle(yTitle.c_str());
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.0);

  h1->SetLineWidth(2);
  h1->SetLineColor(kBlue);

}

/*
void HistogramBuilder::SetAxises(TH1F * h1, string xTitle, string yTitle){
  // X Axis
  h1->GetXaxis()->SetTitle(xTitle.c_str());
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitleOffset(0.9);
  
  // Y Axis
  h1->GetYaxis()->SetTitle(yTitle.c_str());
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.0);

  h1->SetLineWidth(2);
  h1->SetLineColor(kBlue);

}
void HistogramBuilder::SetAxises(TH2F * h2, string xTitle, string yTitle){
  // X Axis
  h2->GetXaxis()->SetTitle(xTitle.c_str());
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.9);

  // Y Axis
  h2->GetYaxis()->SetTitle(yTitle.c_str());
  h2->GetYaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleOffset(1.0);

  h2->SetLineWidth(2);
  h2->SetLineColor(kBlue);
}
*/


/*
 * Iterates through all the stored histograms and writes them
 * to a root file.
 */

void HistogramBuilder::Write(){
  for( map<string,TH1 *>::iterator it = histograms_.begin(); it!= histograms_.end(); ++it){
    it->second->Write();
  }
}


/*
 * Handles wrapping between two angles.
 * So they are never more than 2 PI radians separated.
 * Returns values between -PI and PI.
 */
/*

float HistogramBuilder::WrapCheck(float phi1, float phi2){
  static const float PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
  
  float delta_phi = phi1 - phi2;
  if (delta_phi > PI) delta_phi-= 2*PI;
  if(delta_phi < -1*PI) delta_phi += 2*PI;
    
  return delta_phi;
};

*/
