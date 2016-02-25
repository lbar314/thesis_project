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
  float alpha; // alpha is the angle in y-radius plane in local frame
  float beta;  // beta is the angle in xz plane, taken from z axis, growing counterclockwise
  int nRowPatt;
  int nColPatt;
} clSumm;

TObjArray arrMCTracks; // array of hit arrays for each particle

void debug(int nev=-1)
{
  TFile* boh = TFile::Open("boh.root","recreate");
  TObjArray bitsArr;
  vector<pair<int,Topology> > mappa;
  vector<pair<int,int> > clashes;
  ofstream a("list.txt");
  ofstream b("prova.txt");
  clSumm cSum;
  int pippobaudo=0;
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

  printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  TopDatabase DB;

  for (int iEvent = 0; iEvent < ntotev; iEvent++) {
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

          /*double PG[3] = {(double)pHit->GetPXG(), (double)pHit->GetPYG(), (double)pHit->GetPZG()}; //Momentum at hit-point in Global Frame
          double PL[3];
          if (TMath::Abs(PG[0])<10e-7 && TMath::Abs(PG[1])<10e-7) {
            pHit->Dump();
            int lb = pHit->GetTrack();
            stack->Particle(lb)->Print();
            continue;
          }
          mat->MasterToLocalVect(PG,PL); //Momentum in local Frame
          //printf(">> %e %e   %e %e   %e %e\n",PG[0],PL[0],PG[1],PL[1],PG[2],PL[2]);*/

          Double_t alpha1 = TMath::ACos(TMath::Abs(dirHit[1])/TMath::Sqrt(dirHit[0]*dirHit[0]+dirHit[1]*dirHit[1]+dirHit[2]*dirHit[2])); //Polar Angle
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
	        DB.AccountTopology(*cl, cSum.dX, cSum.dZ, cSum.alpha, cSum.beta);

          /*
            Topology top(*cl);
            int hash = top.GetHash();
            if(primo<6){
              cout << "*****************************\n";
              Topology::printCluster(*cl,cout);
              cout << "Checking print\n";
              top.printTop(cout);
              printf("xCOG: %f + (%f) xCOG: %f + (%f) fired: %d\n", top.GetxCOGPix(),top.GetxCOGshift(),top.GetzCOGPix(),top.GetzCOGshift(), top.GetFiredPixels());
              primo++;
            }
          */
          vector<int> v;
          for(int i=0; i<10; i++){
            Topology top(*cl);
            if(top.GetHash()==-2142478390){
              cout << "Something wrong" << endl;
              exit(1);
            }
            v.push_back(Topology::FuncMurmurHash2(top.GetPattern().data(),(int)top.GetPattern().length()));
          }
          bool er = false;
          for(int j=0; j<v.size(); j++){
            for(int k=j+1; k<v.size(); k++){
              if(v[j]!=v[k]) er = true;
            }
          }
          if(er) for(int d=0; d<v.size(); d++) cout << v[d] << endl;
          //
          //a << ilr << " " << modID << " " << col << " " << row << " " << hash << endl;
          /*
            bool newTop = true;
            for (int i = 0; i < mappa.size(); ++i) {
              if (mappa[i].first == hash) {
                newTop = false;
                if (mappa[i].second.GetPattern() != top.GetPattern() || mappa[i].second.GetUniqueID() != top.GetUniqueID()) {
                  bool newClash = true;
                  for (int j = 0; j < clashes.size(); ++j) {
                    if(clashes[j].first == hash) {
                      clashes[j].second++;
                      newClash = false;
                      break;
                    }
                  }
                  if (newClash) clashes.push_back(pair<int,int>(hash,1));
                }
              }
            }
            if (newTop) mappa.push_back(pair<int,Topology>(hash,top));
          */
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
          /*
          if (clsize==5) {
            printf("\nL%d(%c) Mod%d, Cl:%d | %+5.1f %+5.1f (%d/%d)|H:%e %e %e | C:%e %e %e\n",ilr,cl->TestBit(kSplit) ? 'S':'N',
             modID,icl,(txyzH[0]-xyzClTr[0])*1e4,(txyzH[2]-xyzClTr[2])*1e4, row,col,
             gxyzH[0],gxyzH[1],gxyzH[2],xyzClGlo[0],xyzClGlo[1],xyzClGlo[2]);
            cl->Print();
            pHit->Print();
            //
            double a0,b0,c0,a1,b1,c1,e0;
            pHit->GetPositionL0(a0,b0,c0,e0);
            pHit->GetPositionL(a1,b1,c1);
            float cloc[3];
            cl->GetLocalXYZ(cloc);
            printf("LocH: %e %e %e | %e %e %e\n",a0,b0,c0,a1,b1,c1);
            printf("LocC: %e %e %e | %e %e %e\n",cloc[0],cloc[1],cloc[2],xyzClTr[0],xyzClTr[1],xyzClTr[2]);
          }
          */
          //
        }
      }
    }

    //    layerClus.Clear();
    //
    arrMCTracks.Delete();
  }//event loop
  cout<<"Number of errors in conversion from topology to word: " << pippobaudo << endl;
  boh->cd();
  boh->WriteObject(&bitsArr,"bitsArr","kSingleKey");
  boh->Close();
  delete boh;
  arrMCTracks.Delete();
  //
  DB.EndAndSort();
  DB.SetThresholdCumulative(0.95);
  cout << "Over threshold: : "<< DB.GetOverThr()<<endl;
  DB.Grouping(10,10);
  DB.BuildMap();
  DB.PrintDB("Database1.txt");
  TFile* flDB = TFile::Open("TopologyDatabase.root", "recreate");
  flDB->WriteObject(&DB,"DB","kSingleKey");
  flDB->Close();
  delete flDB;

  a.close();
  b.close();
}
