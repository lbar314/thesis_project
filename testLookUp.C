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

typedef struct {
  int evID;
  int volID;
  int lrID;
  int clID;
  int nPix;
  int nX;
  int nZ;
  int q;
  float pt;
  float eta;
  float phi;
  float xyz[3];
  float dX;
  float dY;
  float dZ;
  bool split;
  bool prim;
  int  pdg;
  int  ntr;
  float alpha; // alpha is the angleTopDat in y-radius plane in local frame
  float beta;  // beta is the angle in xz plane, taken from z axis, growing counterclockwise
  int nRowPatt;
  int nColPatt;
} clSumm;

TObjArray arrMCTracks; //array of hit arrays for each particle

void testLookUp(string inputfile="dizionario.txt", int repetitions = 30, int nev=-1, string outputfile="timeLookUp.txt"){

  LookUp finder(inputfile);
  ofstream time_output(outputfile);
  TStopwatch timerLookUp;
  //
  int num = finder.GetOver() + 49;
  TH1F* hCheck = new TH1F("hCheck","hCheck",num,-0.5,num-0.5);
  hCheck->SetDirectory(0);
  hCheck->SetFillColor(kBlue);
  hCheck->SetFillStyle(3005);

  clSumm cSum;
  int primo=0;

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

  Printf("N Events : %i \n",ntotev+1);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  for(int iRep=0; iRep<repetitions; iRep++){
    for (int iEvent = 0; iEvent < ntotev; iEvent++) {
      Printf("\n Rep %i / %i Event %i / %i\n",iRep+1,repetitions,iEvent+1,ntotev);
      Int_t totClusters=0;
      runLoader->GetEvent(iEvent);
      AliStack *stack = runLoader->Stack();
      cluTree=dl->TreeR();
      hitTree=dl->TreeH();
      hitTree->SetBranchAddress("ITS",&hitList);
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
      for(int iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop of the hits
        hitTree->GetEntry(iEnt);
        int nh = hitList->GetEntries();
        for(int iHit=0; iHit<nh;iHit++){
          AliITSMFTHit *pHit = (AliITSMFTHit*)hitList->At(iHit);
          int mcID = pHit->GetTrack();
  	      //Printf("MCid: %d %d %d Ch %d\n",iEnt,iHit, mcID, pHit->GetChip());
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
      //Printf(" tree entries: %lld\n",cluTree->GetEntries());
      //
      for (int ilr=0;ilr<nlr;ilr++) {
        AliITSURecoLayer* lr = its->GetLayerActive(ilr);
        TClonesArray* clr = lr->GetClusters();
        int nClu = clr->GetEntries();
        //Printf("Layer %d : %d clusters\n",ilr,nClu);
        //
        for (int icl=0;icl<nClu;icl++) {
          if(icl%100==0)printf("ilr: %d icl: %d / %d\r", ilr, icl, nClu);
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
            //  Printf("To check for %d (mod:%d) N=%d from %d\n",icl,modID,nclSn,offs);
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
                    Printf("Discard clusters of module %d:\n",modID);
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
          if (!mat) {Printf("failed to get matrix for module %d\n",cl->GetVolumeId());}
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
  	        //Printf("check %d/%d LB %d  %p\n",il,nLab,labels[il],htArr);
  	        if (!htArr) {Printf("did not find MChits for label %d ",labels[il]); cl->Print(); continue;}
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
              Printf("did not find MChit for label %d on module %d ",il,modID);
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

            Double_t alpha1 = TMath::ACos(TMath::Abs(dirHit[1])/TMath::Sqrt(dirHit[0]*dirHit[0]+dirHit[1]*dirHit[1]+dirHit[2]*dirHit[2]));          //Polar Angle
            float alpha2 = (float) alpha1; //convert to float
            cSum.alpha = alpha2;

            Double_t beta1;
            beta1 = TMath::ATan2(dirHit[0],dirHit[2]); //Azimuthal angle, values from -Pi to Pi
            float beta2 = (float) beta1;
            cSum.beta = beta2;

            cSum.evID = iEvent;
            cSum.volID = cl->GetVolumeId();
            cSum.lrID = ilr;
            cSum.clID = icl;
            cSum.nPix = cl->GetNPix();
            cSum.nX   = cl->GetNx();
            cSum.nZ   = cl->GetNz();
            cSum.q    = cl->GetQ();
            cSum.split = cl->TestBit(kSplit);
            cSum.dX = (txyzH[0]-xyzClTr[0])*1e4;
            cSum.dY = (txyzH[1]-xyzClTr[1])*1e4;
            cSum.dZ = (txyzH[2]-xyzClTr[2])*1e4;
            cSum.nRowPatt = cl-> GetPatternRowSpan();
            cSum.nColPatt = cl-> GetPatternColSpan();

            int a;
            bool restart = (totClusters==0 && iEvent==0) ? true : false;
            timerLookUp.Start(restart);
            a=finder.GroupFinder(*cl);
            timerLookUp.Stop();
            totClusters++;
            if(iRep==0) hCheck->Fill(a);
            //
            int label = cl->GetLabel(0);
            TParticle* part = 0;
            if (label>=0 && (part=stack->Particle(label)) ) {
              cSum.pdg = part->GetPdgCode();
              cSum.eta = part->Eta();
              cSum.pt  = part->Pt();
              cSum.phi = part->Phi();
              cSum.prim = stack->IsPhysicalPrimary(label);
            }
            cSum.ntr = 0;
            for (int ilb=0;ilb<3;ilb++) if (cl->GetLabel(ilb)>=0) cSum.ntr++;
            for (int i=0;i<3;i++) cSum.xyz[i] = xyzClGloF[i];
            //
          }
        }
      }
      //    layerClus.Clear();
      //
      arrMCTracks.Delete();
    }//event loop
    time_output << timerLookUp.RealTime()/ntotev << " " << timerLookUp.CpuTime()/ntotev << endl;
  }
  arrMCTracks.Delete();
  //
  time_output.close();
  cout << endl;
  //
  TCanvas* c = new TCanvas("c","c");
  c->SetLogx();
  c->SetLogy();
  hCheck->Scale(1./hCheck->Integral());
  hCheck->Draw();
  TFile* ciccio = TFile::Open("./histos.root","UPDATE");
  ciccio->WriteObject(hCheck,"lookup","kSingleKey");
  ciccio->Close();
}
