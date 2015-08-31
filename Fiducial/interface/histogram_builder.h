#ifndef __HISTOGRAMBUILDER_H__
#define __HISTOGRAMBUILDER_H__

/*
 * Author Chris Anelli
 * 8.20.2015
 */

#include "TH1F.h"
#include "TH2.h"
#include <map>

using namespace std;


class HistogramBuilder{
  
  //int mynumber = 13;
  
 public:

  //HistogramBuilderTwo(); //Constructor
     
  /*
   * Functions for HistogramBuilderTwo                                       
   */
  
  void FillCutFlowHistograms(string prefix, int cut_step, double weight=1);
  void fillPtHistograms(string prefix, float pt,  double weight=1);
  void fillEtaHistograms(string prefix, float eta, double weight =1);
  
  //void fillEtaPhiHistograms(float eta, float phi, string key, double weight=1);
  //void fillDeltaEtaDeltaPhiHistograms(float eta1, float eta2, 
  //float phi1, float phi2,  string key, double weight=1);
  //void fillCountHistogram(string key,double weight=1);    
  //void fillWeightHistograms(float weight_val, string key, double weight=1);
  //void fillScatterPt(float pt1, float pt2, string key, double weight);
  
  map<string, TH1 *> GetHistograms() { return _histograms;};
  void Write();
  
  //float WrapCheck(float phi1, float phi2);
  //int getMyNumber(); //{return mynumber;}

 private:
  
  map<string, TH1*> _histograms;

  //map<string, TH1F*> _h1CutFlow;
  map<string,TH1F*> _h1Pt;
  map<string,TH1F*> _h1InvPt;
  map<string,TH2F*> _h2ScatterPt;
  map<string,TH1F*> _h1Energy;                                 
  map<string,TH1F*> _h1Eta;
  map<string,TH1F*> _h1Phi;
  map<string,TH2F*> _h2EtaPhiMap;
  map<string,TH1F*> _h1DeltaEta;                                   
  map<string,TH1F*> _h1DeltaPhi;
  map<string,TH2F*> _h2DeltaEtaDeltaPhi;
  map<string,TH1F*> _h1Trig;
  map<string,TH1F*> _h1Weight;
  map<string,TH1F*> _h1Counter; 
  
  /*
   * Helper Functions
   */

  void SetAxises(TH1 * h1, string xtitle, string ytitle);

  //void SetAxises(TH1F * h1, string xtitle, string ytitle);
  //void SetAxises(TH2F * h2, string xtitle, string ytitle);
  
};

#endif
