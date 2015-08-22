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

void HistogramBuilder::FillCutFlowHistograms(string key, int cut_step, double weight){
  if(!_histograms.count(key)){
    _histograms[key] = new TH1F(Form("%s_CutFlow",key.c_str()),
			       Form("%s CutFlow",key.c_str()),
			       10, 0, 10);
    _histograms[key]->GetYaxis()->SetTitle("CutFlow");
    _histograms[key]->Sumw2();
  }
  _histograms[key]->Fill(cut_step, weight);
}

/*
 *Pt Histograms
 */
void HistogramBuilder::fillPtHistograms(string key, float pt, double weight){
  if(!_histograms.count(key)){
        
    _histograms[key] = new TH1F(Form("%s_Pt",key.c_str()),
					  Form("%s Pt",key.c_str()),
					  100,0,200);
    SetAxises(_histograms[key],"Pt (GeV)", "Counts");
    _histograms[key]->Sumw2();
  } 
  _histograms[key]->Fill(pt, weight);
}


/*
 *Counting Histograms
 *Fills the 1 bin.
 */
void HistogramBuilder::fillCountHistogram(string key, double weight){
  if(!_h1Counter.count(key)){                                                   
    _h1Counter[key] = new TH1F(Form("%s_Count",key.c_str()),    
                                               Form("%s Count",key.c_str()),    
                                               2, 0, 2);
    _h1Counter[key]->GetYaxis()->SetTitle("Counts");
    _h1Counter[key]->Sumw2();
  }                                                                             
  _h1Counter[key]->Fill(1, weight);
}
                                                                               

/*      
 *Trigger Histograms
 */
void HistogramBuilder::fillTrigHistograms(bool trigDecision,string key, double weight){ 
  if(!_h1Trig.count(key)){                                                   
    _h1Trig[key] = new TH1F(Form("%s_Trig",key.c_str()),     
                                            Form("%s Trigger",key.c_str()),  
                                            2, 0, 2);
    SetAxises(_h1Trig[key],"Trigger Decision", "Counts");
    _h1Trig[key]->Sumw2();
  }                                                      
                       
  _h1Trig[key]->Fill(trigDecision, weight);                                        
     
}

/*      
 *Weight Histograms
 */
void HistogramBuilder::fillWeightHistograms(float weight_val,string key, double weight){ 
  if(!_h1Weight.count(key)){                                                   
    _h1Weight[key] = new TH1F(Form("%s_Weight",key.c_str()),     
					      Form("%s Weight",key.c_str()),  
					      15000, 0, 150);
    SetAxises(_h1Weight[key],"Weight", "Counts");
  }                                                      
                       
  _h1Weight[key]->Fill(weight_val, weight);                                        
     
}
      

/*                                          
 *Energy Histograms
 */

void HistogramBuilder::fillEnergyHistograms(float energy, string key, double weight){  
  if(!_h1Energy.count(key)){                                             
    _h1Energy[key] = new TH1F(Form("%s_Energy",key.c_str()),
                                            Form("%s Energy",key.c_str()), 
                                            2100, -5.0, 100.0);
    SetAxises(_h1Energy[key],"Energy (GeV)", "Counts");
    _h1Energy[key]->Sumw2();
  }                                                                          
  _h1Energy[key]->Fill(energy, weight);                                            
}  

/*                                                                             
 *Eta Phi Histograms                                                         
 */

void HistogramBuilder::fillEtaPhiHistograms(float eta, float phi, string key, double weight){
  if(!_h1Eta.count(key)){
    _h1Eta[key] = new TH1F(Form("%s_Eta",key.c_str()),
                                           Form("%s Eta",key.c_str()),
                                           500, -1.5, 1.5);  //HO has 72 iphis and 30 ietas
    SetAxises(_h1Eta[key], "#eta", "Counts");
    _h1Eta[key]->Sumw2();
    
  }
  _h1Eta[key]->Fill(eta, weight);

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

/*
 *Delta Eta Delta Phi Histograms
 */

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


/*
if(!_h1InvPt.count(key)){
        
    _h1InvPt[key] = new TH1F(Form("%s_InvPt",key.c_str()),
					  Form("%s Inv Pt",key.c_str()),
					  800,-1,1);
    SetAxises(_h1InvPt[key], "#frac{1}{P_t}", "Counts");
    _h1InvPt[key]->Sumw2();
  } 
 _h1InvPt[key]->Fill(1.0/pt, weight);
*/

/*
 * Scatter Pt Plot
 */

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
 * Handles wrapping between two angles.
 * So they are never more than 2 PI radians separated.
 * Returns values between -PI and PI.
 */



float HistogramBuilder::WrapCheck(float phi1, float phi2){
  static const float PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
  
  float delta_phi = phi1 - phi2;
  if (delta_phi > PI) delta_phi-= 2*PI;
  if(delta_phi < -1*PI) delta_phi += 2*PI;
    
  return delta_phi;
  /*
  float delta_case;
   delta_case = delta_phi;
  if(delta_phi < -PI){
    delta_case = (2*PI + delta_phi);
  }
  if(delta_phi > PI){
    delta_case = (delta_phi - 2*PI);
  }
    
  //float delta_mod = fmod(pos_angle + PI, 2*PI) - PI;
 
  // Use the float modulus operator
  */
};
