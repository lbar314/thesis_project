#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "AliITSMFTClusterPix.h"
#include "AliITSURecoLayer.h"
#include "AliITSURecoDet.h"
#include "AliITSMFTHit.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliGeomManager.h"
#include "AliStack.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGeoMatrix.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include <Riostream.h>
#include "./MinimTopology.h"
#include "./Dictionary.h"
#include "BuildDictionary.h"
#include "LookUp.h"
#include <map>
#include <math.h>

#endif

using namespace std;

enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

void timeDistr(string inputfile="dizionario.txt", int repetitions = 1, int nev=-1){

  LookUp finder(inputfile);

  TStopwatch timerLookUp;

  TH1F* hr = new TH1F("hr","LookUp: real time",100,0,1.5e-6);
  hr->SetDirectory(0);
  hr->GetXaxis()->SetTitle("t (s)");
  hr->SetFillColor(kRed);
  hr->SetFillStyle(3005);

  TH1F* hc = new TH1F("hc","LookUp: CPU time",100,0,1.5e-6);
  hc->SetDirectory(0);
  hc->GetXaxis()->SetTitle("t (s)");
  hc->SetFillColor(kRed);
  hc->SetFillStyle(3005);

  const int kSplit=0x1<<22;
  const int kSplCheck=0x1<<23;
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Align/Data",Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Calib/RecoParam",Form("local://%s",gSystem->pwd()));
  man->SetRun(0);

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();
  runLoader->LoadSDigits();
  runLoader->LoadHits();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  TObjArray algITS;
  AliGeomManager::LoadAlignObjsFromCDBSingleDet("ITS",algITS);
  AliGeomManager::ApplyAlignObjsToGeom(algITS);
  //
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(true);
  AliITSMFTClusterPix::SetGeom(gm);
  //
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
  //
  TTree *cluTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSMFTHit");
  //
  int nlr=its->GetNLayersActive();
  int ntotev = (int)runLoader->GetNumberOfEvents();

  Printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  for(int iRep=0; iRep<repetitions; iRep++){
    for (int iEvent = 0; iEvent < ntotev; iEvent++){
      Printf("\n Rep %i / %i Event %i / %i\n",iRep+1,repetitions,iEvent+1,ntotev);
      runLoader->GetEvent(iEvent);
      AliStack *stack = runLoader->Stack();
      cluTree=dl->TreeR();
      //
      // read clusters
      for (int ilr=nlr;ilr--;) {
        TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",ilr));
        if (!br) {Printf("Did not find cluster branch for lr %d\n",ilr); exit(1);}
        br->SetAddress(its->GetLayerActive(ilr)->GetClustersAddress());
      }
      cluTree->GetEntry(0);
      its->ProcessClusters();
      //
      // read hits
      //
      for (int ilr=0;ilr<nlr;ilr++) {
        AliITSURecoLayer* lr = its->GetLayerActive(ilr);
        TClonesArray* clr = lr->GetClusters();
        int nClu = clr->GetEntries();
        //Printf("Layer %d : %d clusters\n",ilr,nClu);
        //
        for (int icl=0;icl<nClu;icl++) {
          AliITSMFTClusterPix *cl = (AliITSMFTClusterPix*)clr->At(icl);
          timerLookUp.Start();
          finder.GroupFinder(*cl);
          timerLookUp.Stop();
          hr->Fill(timerLookUp.RealTime());
          hc->Fill(timerLookUp.CpuTime());
        }
      }
    }//event loop
  }
  //
  cout << endl;
  TCanvas* c = new TCanvas("c","c");
  c->Divide(2,1);
  c->cd(1);
  hr->Draw();
  c->cd(2);
  hc->Draw();
}
