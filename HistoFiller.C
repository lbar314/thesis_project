#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#endif

using std::string;
using std::ifstream;

void HistoFiller(string filename){
  TH1F* hReal = new TH1F("hReal", "Look-up: Real time",100,-1e-7,2e-1);
  hReal->SetFillColor(kRed);
  hReal->SetFillStyle(3005);
  hReal->GetXaxis()->SetTitle("t(s)");
  //
  TH1F* hCpu = new TH1F("hCpu", "Look-up: CPU time",100,-1e-7,2e-1);
  hCpu->SetFillColor(kBlue);
  hCpu->SetFillStyle(3005);
  hCpu->GetXaxis()->SetTitle("t(s)");
  //
  double cpuvalue;
  double realvalue;
  //
  ifstream in(filename);
  if(!in.is_open()){
    std::cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in >> realvalue >> cpuvalue){
      hCpu->Fill(cpuvalue);
      hReal->Fill(realvalue);
    }
  }
  in.close();

  TCanvas* c = new TCanvas("c","");
  c->Divide(2,1);
  c->cd(1);
  hReal->Draw();
  c->cd(2);
  hCpu->Draw();
}
