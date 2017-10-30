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
#include "./cluster2string.h"

#endif

using namespace std;

enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

void StoreClusters(int nev=-1, string outputfile="clusterlist.bin"){

  ofstream check("store.txt");

  ofstream file_output(outputfile, ios::out | ios::binary);
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

  int ntotclusters=0;

  Printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  for (int iEvent = 0; iEvent < ntotev; iEvent++) {
    Printf("Event %i / %i\n",iEvent+1,ntotev);
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
    for (int ilr=0;ilr<nlr;ilr++) {
      AliITSURecoLayer* lr = its->GetLayerActive(ilr);
      TClonesArray* clr = lr->GetClusters();
      int nClu = clr->GetEntries();
      //Printf("Layer %d : %d clusters\n",ilr,nClu);
      //
      for (int icl=0;icl<nClu;icl++){
        AliITSMFTClusterPix *cl = (AliITSMFTClusterPix*)clr->At(icl);
        string str;
        FromCluster2String(*cl,str);
        int stringSize = (int)(str.size());
        check << ntotclusters << ") string size: " << stringSize << endl;
        MinimTopology a(str);
        check << a;
        file_output.write(reinterpret_cast<char *>(&stringSize), sizeof(int));
        file_output.write(str.c_str(), stringSize);
        ntotclusters++;
      }
    }
  }//event loop
  //
  check.close();
  cout << "Number of clusters: " << ntotclusters << endl;
  file_output.close();
  cout << endl;
}
