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
#include "./LookUp.h"
#include <map>

#endif

using namespace std;

enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

TObjArray arrMCTracks; // array of hit arrays for each particle

void testLookUp(string inputfile,int nRepetintions=100,int nev=-1){

  LookUp finder(inputfile);
  TStopwatch timerLookUp;
  //
  TCanvas* c = new TCanvas("c","cTime");
  c->Divide(2,1);
  //
  TH1F* timerReal = new TH1F("timerReal","Real time with the old method",50,0,1e-1);
  timerReal->SetDirectory(0);
  timerReal->GetXaxis()->SetTitle("t (s)");
  timerReal->SetFillColor(kBlue);
  timerReal->SetFillStyle(3008);
  timerReal->SetNdivisions(505,"X");
  //
  TH1F* timerCpu = new TH1F("timerCpu","CPU time with the old method",50,0,1e-1);
  timerCpu->SetDirectory(0);
  timerCpu->GetXaxis()->SetTitle("t (s)");
  timerCpu->SetFillColor(kRed);
  timerCpu->SetFillStyle(3008);

  //
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
  Double_t xg1,yg1,zg1=0.,xg0,yg0,zg0=0.,tg0;
  Double_t xExit,yExit,zExit,xEnt,yEnt,zEnt,tof1;

  //
  TTree *cluTree = 0x0;
  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSMFTHit");
  //
  float xyzClGloF[3];
  Double_t xyzClGlo[3],xyzClTr[3];
  int labels[3];
  int nLab = 0;
  int nlr=its->GetNLayersActive();
  int ntotev = (int)runLoader->GetNumberOfEvents();

  printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  for(int iRep=0; iRep<nRepetintions; iRep++){
    for (int iEvent = 0; iEvent < ntotev; iEvent++) {
      Int_t totClusters=0;
      printf("\n Event %i \n",iEvent);
      runLoader->GetEvent(iEvent);
      AliStack *stack = runLoader->Stack();
      cluTree=dl->TreeR();
      hitTree=dl->TreeH();
      hitTree->SetBranchAddress("ITS",&hitList);
      //
      // read clusters
      for (int ilr=nlr;ilr--;) {
        TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",ilr));
        if (!br) {printf("Did not find cluster branch for lr %d\n",ilr); exit(1);}
        br->SetAddress(its->GetLayerActive(ilr)->GetClustersAddress());
      }
      cluTree->GetEntry(0);
      its->ProcessClusters();
      //
      // read hits
      for(int iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop of the hits
        hitTree->GetEntry(iEnt);
        int nh = hitList->GetEntries();
        for(int iHit=0; iHit<nh;iHit++){
          AliITSMFTHit *pHit = (AliITSMFTHit*)hitList->At(iHit);
          int mcID = pHit->GetTrack();
  	      //printf("MCid: %d %d %d Ch %d\n",iEnt,iHit, mcID, pHit->GetChip());
          TClonesArray* harr = arrMCTracks.GetEntriesFast()>mcID ? (TClonesArray*)arrMCTracks.At(mcID) : 0;
          if (!harr) {
            harr = new TClonesArray("AliITSMFTHit"); // 1st encounter of the MC track
            arrMCTracks.AddAtAndExpand(harr,mcID);
          }
          //
          new ( (*harr)[harr->GetEntriesFast()] ) AliITSMFTHit(*pHit);
        }
      }
      //    return;
      //
      // compare clusters and hits
      //
      printf(" tree entries: %lld\n",cluTree->GetEntries());
      //
      for (int ilr=0;ilr<nlr;ilr++) {
        AliITSURecoLayer* lr = its->GetLayerActive(ilr);
        TClonesArray* clr = lr->GetClusters();
        int nClu = clr->GetEntries();
        //printf("Layer %d : %d clusters\n",ilr,nClu);
        //
        for (int icl=0;icl<nClu;icl++) {
          printf("iRep: %d iEv: %d ilr: %d icl: %d / %d\n", iRep, iEvent, ilr, icl, nClu);
          AliITSMFTClusterPix *cl = (AliITSMFTClusterPix*)clr->At(icl);
          int modID = cl->GetVolumeId();

          //------------ check if this is a split cluster
          int sInL = modID - gm->GetFirstChipIndex(ilr);
          if (!cl->TestBit(kSplCheck)) {
            cl->SetBit(kSplCheck);
            // check if there is no other cluster with same label on this module
            AliITSURecoSens* sens = lr->GetSensor(sInL);
            int nclSn = sens->GetNClusters();
            int offs = sens->GetFirstClusterId();
            //  printf("To check for %d (mod:%d) N=%d from %d\n",icl,modID,nclSn,offs);
            for (int ics=0;ics<nclSn;ics++) {
              AliITSMFTClusterPix* clusT = (AliITSMFTClusterPix*)lr->GetCluster(offs+ics); // access to clusters
              if (clusT==cl) continue;
              for (int ilb0=0;ilb0<3;ilb0++) {
                int lb0 = cl->GetLabel(ilb0); if (lb0<=-1) break;
                for (int ilb1=0;ilb1<3;ilb1++) {
                  int lb1 = clusT->GetLabel(ilb1); if (lb1<=-1) break;
                  if (lb1==lb0) {
                    cl->SetBit(kSplit);
                    clusT->SetBit(kSplit);
                    //printf("Discard clusters of module %d:\n",modID);
                    //cl->Print();
                    //clusT->Print();
                    break;
                  }
                }
              }
            }
          }
          timerLookUp.Start(!totClusters);
          finder.GroupFinder(*cl);
          timerLookUp.Stop();
          totClusters++;
        }
      }
      //    layerClus.Clear();
      //
      arrMCTracks.Delete();
      timerCpu->Fill(timerLookUp.CpuTime());
      timerReal->Fill(timerLookUp.RealTime());
    }//event loop
    arrMCTracks.Delete();
  }//repetition loop
  c->cd(1);
  timerReal->Draw();
  c->cd(2);
  timerCpu->Draw();
  c->Print("Time.pdf");
}
