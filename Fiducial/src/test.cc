/*
 * Quick Test File to see if my makefile
 * is able to access the root libraries
 */
#include "TFile.h"
#include "TH1F.h"
#include <iostream>

using namespace std;

int main(){
  TFile * outfile = new TFile("test.root", "RECREATE");
  cout << "hello world" << endl;
  TH1F * h1 = new TH1F("h1", "h1", 10, 0, 10);
  
  h1->Write();
  outfile->Close();

}
