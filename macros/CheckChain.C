#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <Riostream.h>
#include <math.h>

#endif

void CheckChain(string inputfile="histos.root"){
  TFile* ciccio = TFile::Open(inputfile.data());
  TH1F* h1 = (TH1F*)ciccio->Get("costr");
  TH1F* h2 = (TH1F*)ciccio->Get("lookup");
  TCanvas* c = new TCanvas("c","c");
  c->Divide(2,1);
  c->cd(1);
  h1->Draw();
  c->cd(2);
  h2->Draw();
  int n1 = h1->GetNbinsX();
  int n2 = h2->GetNbinsX();
  if((n1-n2)!=0){
    cout << "different number of GroupIDs" << endl;
    return;
  }
  bool same = true;
  for(int i =0; i<=n1+1; i++){
    if(fabs(h1->GetBinContent(i)-h2->GetBinContent(i))>(h1->GetBinContent(i))*1.e-3){
      same = false;
      break;
    }
  }
  if(same) cout <<"Same distributions! :D" << endl;
  else cout <<"Different distributions :(" << endl;
}
