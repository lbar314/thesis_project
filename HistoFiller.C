#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#endif

void HistoFiller(string inputfile="timeLookUp.txt"){

  TH1F* hr = new TH1F("hr","LookUp: real time",100,0,2e-1);
  hr->SetDirectory(0);
  hr->GetXaxis()->SetTitle("t (s)");
  hr->SetFillColor(kRed);
  hr->SetFillStyle(3005);

  TH1F* hc = new TH1F("hc","LookUp: CPU time",100,0,2e-1);
  hc->SetDirectory(0);
  hc->GetXaxis()->SetTitle("t (s)");
  hc->SetFillColor(kRed);
  hc->SetFillStyle(3005);

  double real;
  double cpu;

  std::ifstream in(inputfile);
  if(!in.is_open()){
    std::cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in >> real >> cpu){
      hr->Fill(real);
      hc->Fill(cpu);
    }
  }
  in.close();

  TCanvas* c = new TCanvas("c","c");
  c->Divide(2,1);
  c->cd(1);
  hr->Draw();
  c->cd(2);
  hc->Draw();
}
