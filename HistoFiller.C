#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#endif

void HistoFiller(string inputfile="timeLookUp.txt"){

  TH1F* hr = new TH1F("hr","LookUp: real time",100,0,1.5e-1);
  hr->SetDirectory(0);
  hr->GetXaxis()->SetTitle("t (s)");
  hr->SetFillColor(kRed);
  hr->SetFillStyle(3005);

  TH1F* hc = new TH1F("hc","LookUp: CPU time",100,0,1.5e-1);
  hc->SetDirectory(0);
  hc->GetXaxis()->SetTitle("t (s)");
  hc->SetFillColor(kRed);
  hc->SetFillStyle(3005);

  double real;
  double cpu;

  double realmax=0.;
  double realmin=0.;
  double cpumax=0.;
  double cpumin=0.;

  double cpumean=0;
  double cpuerr=0;
  double realmean=0;
  double realerr=0;

  int count = 0;

  std::ifstream in(inputfile);
  if(!in.is_open()){
    std::cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in >> real >> cpu){
      count++;
      if(real>realmax) realmax = real;
      if(count==1) realmin = real;
      else {
        if(real<realmin) realmin = real;
      }
      realmean += real;
      if(cpu>cpumax) cpumax = cpu;
      if(count==1) cpumin = cpu;
      else {
        if(cpu<cpumin) cpumin = cpu;
      }
      cpumean += cpu;
      hr->Fill(real);
      hc->Fill(cpu);
    }
  }
  in.close();
  cpumean /= count;
  realmean /= count;
  realerr = (realmax-realmin)/2;
  cpuerr = (cpumax-cpumin)/2;

  cout << "real: " << realmean << " +/- " << realerr << endl;
  cout << "cpu: " << cpumean << " +/- " << cpuerr << endl;

  TCanvas* c = new TCanvas("c","c");
  c->Divide(2,1);
  c->cd(1);
  hr->Draw();
  c->cd(2);
  hc->Draw();
}
