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
#include <Riostream.h>
#include "./Topology.h"
#include "./TopDatabase.h"
#include <map>

#endif

using namespace std;

enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

typedef struct {
  Int_t evID;
  Int_t volID;
  Int_t lrID;
  Int_t clID;
  Int_t nPix;
  Int_t nX;
  Int_t nZ;
  Int_t q;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t xyz[3];
  Float_t dX;
  Float_t dY;
  Float_t dZ;
  Bool_t split;
  Bool_t prim;
  Int_t  pdg;
  Int_t  ntr;
  Float_t alpha; // alpha is the angle in y-radius plane in local frame
  Float_t beta;  // beta is the angle in xz plane, taken from z axis, growing counterclockwise
  Int_t nRowPatt;
  Int_t nColPatt;
} clSumm;

TObjArray arrMCTracks; // array of hit arrays for each particle

void testDecoding5(int nev=-1)
{
  vector<pair<int,Topology> > mappa;
  vector<pair<int,int> > clashes;
  clSumm cSum;

  TopDatabase DB;
  TFile* fl = TFile::Open("TopologyDatabase.root");
  DB = *((TopDatabase*)fl->Get("DB"));
  fl->Close();
  delete fl;
  DB.PrintDB("check.txt");

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
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
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
  Float_t xyzClGloF[3];
  Double_t xyzClGlo[3],xyzClTr[3];
  Int_t labels[3];
  int nLab = 0;
  int nlr=its->GetNLayersActive();
  int ntotev = (Int_t)runLoader->GetNumberOfEvents();

  printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //

  for (Int_t iEvent = 0; iEvent < ntotev; iEvent++) {
    if(iEvent>0) break;
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
    for(Int_t iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop of the hits
      hitTree->GetEntry(iEnt);
      int nh = hitList->GetEntries();
      for(Int_t iHit=0; iHit<nh;iHit++){
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
        if(icl%100 == 0) printf("ilr: %d icl: %d / %d\n", ilr, icl, nClu);
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
                  /*
                  printf("Discard clusters of module %d:\n",modID);
                  cl->Print();
                  clusT->Print();
                  */
                  break;
                }
              }
            }
          }
        }
        //------------
        const AliITSMFTSegmentationPix* segm = gm->GetSegmentation(ilr);
        //
        cl->GetGlobalXYZ(xyzClGloF);
        int clsize = cl->GetNPix();
        for (int i=3;i--;) xyzClGlo[i] = xyzClGloF[i];
        const TGeoHMatrix* mat = gm->GetMatrixSens(modID);
        if (!mat) {printf("failed to get matrix for module %d\n",cl->GetVolumeId());}
        mat->MasterToLocal(xyzClGlo,xyzClTr);
        //
        int col,row;
        segm->LocalToDet(xyzClTr[0],xyzClTr[2],row,col); // effective col/row
        nLab = 0;
        for (int il=0;il<3;il++) {
          if (cl->GetLabel(il)>=0) labels[nLab++] = cl->GetLabel(il);
          else break;
        }
        // find hit info
        for (int il=0;il<nLab;il++) {
          TClonesArray* htArr = (TClonesArray*)arrMCTracks.At(labels[il]);
	        //printf("check %d/%d LB %d  %p\n",il,nLab,labels[il],htArr);
	        if (!htArr) {printf("did not find MChits for label %d ",labels[il]); cl->Print(); continue;}
          //
          int nh = htArr->GetEntriesFast();
          AliITSMFTHit *pHit=0;
          for (int ih=nh;ih--;) {
            AliITSMFTHit* tHit = (AliITSMFTHit*)htArr->At(ih);
            if (tHit->GetChip()!=modID) continue;
            pHit = tHit;
            break;
          }
          if (!pHit) {
            printf("did not find MChit for label %d on module %d ",il,modID);
            cl->Print();
            htArr->Print();
            continue;
          }
          //
          pHit->GetPositionG(xg1,yg1,zg1);
          pHit->GetPositionG0(xg0,yg0,zg0,tg0);
          //
          double txyzH[3],gxyzH[3] = { (xg1+xg0)/2, (yg1+yg0)/2, (zg1+zg0)/2 };
          mat->MasterToLocal(gxyzH,txyzH);

          double rcl = TMath::Sqrt(xyzClTr[0]*xyzClTr[0]+xyzClTr[1]*xyzClTr[1]);
          double rht = TMath::Sqrt(txyzH[0]*txyzH[0]+txyzH[1]*txyzH[1]);
          //
          //Angles determination

          pHit->GetPositionL(xExit,yExit,zExit,gm);
          pHit->GetPositionL0(xEnt,yEnt,zEnt,tof1,gm);

          Double_t dirHit[3]={(xExit-xEnt),(yExit-yEnt),(zExit-zEnt)};

          Topology top(*cl);
          if(top.GetHash()==-2142478390){
            cout << "Something wrong" << endl;
            exit(1);
          }
      	  Int_t num = DB.FromCluster2GroupID(*cl);
        }
      }
    }

    //    layerClus.Clear();
    //
    arrMCTracks.Delete();
  }//event loop
  arrMCTracks.Delete();
}
