#if !defined(__CINT__) || defined(__MAKECINT__)
#include "./Dictionary.h"
#include "./FastSimulation.h"
#include "TH1F.h"
#include "TCanvas.h"
#endif

void testSimulan(unsigned int rep=10000){
  FastSimulation test = FastSimulation("dizionario.txt");
  TH1F* histo = new TH1F("histo","histo",100,-0.5,99.5);
  histo->GetXaxis()->SetTitle("GroupID");
  histo->SetFillStyle(3008);
  histo->SetFillColor(kGreen-3);
  for(unsigned int i=0; i<rep; i++){
    histo->Fill(test.GetRandom());
  }
  TCanvas* c = new TCanvas("c","c");
  c->cd();
  histo->Draw();
}
