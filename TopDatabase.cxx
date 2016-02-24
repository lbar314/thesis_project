#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "./Topology.h"
#include "TFile.h"
#include "./TopDatabase.h"
#include "AliITSMFTClusterPix.h"

ClassImp(TopDatabase)

TopDatabase::TopDatabase():TObject(),fN(0),fTotClusters(0),fThreshold(0),
  fOverThr(0),fNGroups(0),fNmax(1e5),fArrTopologies(),fArrHisto(){
}

TopDatabase::TopDatabase(TopDatabase &ogg):TObject(),fArrHisto(){
  fN=ogg.GetN();
  fArrTopologies.Expand(fN);
  for(int i=0; i<fN; i++){
    fArrTopologies.AddAt(new Topology(*(Topology*)(ogg.GetArrTopologies()->At(i))),i);
  }
  /*
    for(int j=0; j<ogg.GetArrHisto()->GetEntries(); j++){
      fArrHisto.AddAt(new Topology(*(Topology*)(ogg.GetArrTopologies()->At(i))),i);
    }
  */
  fTotClusters = ogg.GetTotClusters();
  fThreshold = ogg.GetThreshold();
  fOverThr = ogg.GetOverThr();
  fNGroups = ogg.GetNGroups();
  fNmax = ogg.GetNmax();
}

TopDatabase::~TopDatabase(){
  fArrTopologies.Delete();
}

void TopDatabase::AccountTopology(const AliITSMFTClusterPix &cluster, Float_t dX, Float_t dZ, Float_t alpha, Float_t beta){

  Topology top(cluster);
  Bool_t newPatt = kTRUE;
  Bool_t Junk = kFALSE;
  fTotClusters++;
  Int_t indTop = -1;
  for(Int_t ip=0;ip<fN;ip++) {
    string pattOld = ((Topology*)fArrTopologies.At(ip))->GetPattern();
    if(top.GetPattern().length() == pattOld.length() && memcmp(top.GetPattern().data(),pattOld.data(),pattOld.length())==0){
      newPatt = kFALSE;
      //cout<< "old" << endl;
      indTop = ip;
      break;
    }
  }
  if(newPatt){
    //cout << "new" << endl;
    if(fN == fNmax){ //Junk bin
      Junk = kTRUE;
    }
    else {
      this->ExpandDB(top);
      indTop=fN-1;
    }
  }
  if(Junk==kTRUE) return;
  ((Topology*)fArrTopologies.At(indTop))->IncreaseCounts();
  ((Topology*)fArrTopologies.At(indTop))->GetHxA().Fill(alpha, dX);
  ((Topology*)fArrTopologies.At(indTop))->GetHzA().Fill(alpha, dZ);
  ((Topology*)fArrTopologies.At(indTop))->GetHxB().Fill(beta, dX);
  ((Topology*)fArrTopologies.At(indTop))->GetHzB().Fill(beta, dZ);
}

void TopDatabase::ExpandDB(const Topology &top){
  fN++;
  fArrTopologies.Expand(fN);
  Topology* top1 = new Topology(top,fN-1);
  fArrTopologies.AddAt(top1,fN-1);
}

void TopDatabase::SetThreshold(Float_t thr){
  fThreshold=thr;
  Int_t nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());
  Int_t over=0;
  Float_t fr=0;
  for(Int_t j=0; j<fN; j++){
    fr=arrFreq[sortIndex[j]];
    if(fr<thr) break;
    else over++;
  }
  fOverThr=over;
}

void TopDatabase::SetThresholdCumulative(Float_t cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  Float_t totFreq = 0.;
  Int_t nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TArrayF provvFreq;
  provvFreq.Set(fN);
  TMath::Sort(fN,arrFreq.GetArray(),sortIndex.GetArray());
  for(Int_t j=0; j<fN; j++){
    provvFreq[j]=arrFreq[sortIndex[j]];
  }
  Int_t over=0;
  while(totFreq < cumulative){
    totFreq+=provvFreq[over++];
  }
  fThreshold=provvFreq[--over];
  while(provvFreq[over]==fThreshold) over++;
  over++;
  fOverThr = over;
}

void TopDatabase::EndAndSort(Int_t mode){
  Topology::SetMode(mode);
  Int_t nPatterns = fN;
  fArrTopologies.Sort();
  Bool_t alreadyProcessed[nPatterns];
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = ((Float_t)(top->GetCounts()))/fTotClusters;
    top->SetFreq(tempFreq);
    arrFreq[i] = tempFreq;
    top->SetFlag(0);
    alreadyProcessed[i]=kFALSE;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());
  for(Int_t i=0; i<nPatterns; i++){
    if(alreadyProcessed[i]==kTRUE) continue;
    else alreadyProcessed[i]=kTRUE;
    Topology* topFreq = (Topology*)fArrTopologies.At(sortIndex[i]);
    topFreq->SetPattID(i);
    Topology* top = (Topology*)fArrTopologies.At(i);
    Int_t refHash = top->GetHash();
    for(Int_t j=i+1; j<nPatterns; j++){
      Topology* top2comp = (Topology*)fArrTopologies.At(j);
      Int_t counter = 0;
      if(top2comp->GetHash()==refHash){
      	alreadyProcessed[j]=kTRUE;
      	if(counter==0){
      	  counter++;
      	  top->SetFlag(1);
      	}
      	top2comp->SetFlag(1);
      }
    }
  }
}

void TopDatabase::PrintDB(const char* output) const{
  ofstream o(output);
  o << "Number of topologies: " << fN << endl << endl;
  for(Int_t i=0; i<fN; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    o << "               POSITION " << i << endl;
    o << Form("\nPattID %4d R:%d C: %d Freq:%f Hash:%d\n\n", top->GetPattID(),top->GetRowSpan(),top->GetColumnSpan(),top->GetFreq(),top->GetHash());
    top->printTop(o);
    o << "\nxCentre: " << top->GetxCOGPix() << " + " << top->GetxCOGshift()
      <<" zCentre: " << top->GetzCOGPix() << " + " << top->GetzCOGshift()
      <<" groupID: " << top->GetGroupID()<< "\n\n...............................................\n";
  }
  o.close();
}

void TopDatabase::Grouping(Int_t NumberofShiftXbins, Int_t NumberofShiftZbins){
  printf("\nGROUPING\n\n");
  Float_t threshold=fThreshold;
  Int_t nPatterns=fN;
  TArrayI shiftXID;
  TArrayI shiftZID;
  Int_t notINgorups=0;
  Float_t shiftXwidth=1./NumberofShiftXbins; //(fraction of the pitch)
  Float_t shiftZwidth=1./NumberofShiftZbins; //(fraction of the pitch)
  Float_t xOffset = shiftXwidth/2;
  Float_t zOffset = shiftZwidth/2;

  shiftXID.Set(nPatterns);
  shiftZID.Set(nPatterns);

  TH1F xShiftDistro("xShiftDistro", "x position", NumberofShiftXbins+1 , -0.5-xOffset, 0.5+xOffset);
  xShiftDistro.SetDirectory(0);

  TH1F zShiftDistro("zShiftDistro", "z position", NumberofShiftZbins+1, -0.5-zOffset, 0.5+zOffset);
  zShiftDistro.SetDirectory(0);

  printf("Assigning shift-IDs:\n");
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    xShiftDistro.Fill(top->GetxCOGshift());
    shiftXID[i]=xShiftDistro.FindBin(top->GetxCOGshift());
    zShiftDistro.Fill(top->GetzCOGshift());
    shiftZID[i]=zShiftDistro.FindBin(top->GetzCOGshift());
  }
  Int_t tempgroupID=0;
  printf("Processing patterns over threshold:\n");
  //first assigning groupID to patterns whose frequency is above the treshold
  for(Int_t t=0; t<nPatterns;t++){
    if(t%100==0)printf("%d / %d\n ", t, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(t);
    if(top->GetFreq()>=threshold){
      top->SetGroupID(tempgroupID);
      tempgroupID++;
      notINgorups++;
    }
  }
  printf("Processing patterns under threshold:\n");
  for(Int_t i=0; i<nPatterns; i++){
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(top->GetGroupID()!=-1) continue;
    top->SetGroupID(tempgroupID);
    //if(i%10000==0)printf("Assigning group ID %d / %d\n ", i, nPatterns);
    for(Int_t j=i+1; j<nPatterns; j++){
      Topology* top1 = (Topology*)fArrTopologies.At(j);
      if(top1->GetGroupID()!=-1) continue;
      else {
      	if(shiftXID[j]==shiftXID[i] && shiftZID[j]==shiftZID[i]
        && top->GetRowSpan()==top1->GetRowSpan() && top->GetColumnSpan()==top1->GetColumnSpan()) top1->SetGroupID(tempgroupID);
      }
    }
    tempgroupID++;
  }
  fNGroups = tempgroupID;
  fOverThr = notINgorups;
  //________________________________OPERATIONS WITH HISTOGRAMS
  fArrHisto.Expand(fNGroups*4);
  fArrHisto.Clear();

  printf("Number of groups: %d\n", fNGroups);

  printf("Summing group-histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    TH2F tempXa(Form("hXA%d",iGroup),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempXb(Form("hXB%d",iGroup),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempZa(Form("hZA%d",iGroup),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempZb(Form("hZB%d",iGroup),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    Bool_t FirstMatch = kTRUE;
    for(Int_t i=0; i<nPatterns; i++){
      // printf("iGroup: %d i: %d\n", iGroup, i);
      Topology* top = (Topology*)fArrTopologies.At(i);
      if(top->GetGroupID()==iGroup){
      	if(FirstMatch==kTRUE){
      	  FirstMatch==kFALSE;
          tempXa = top->GetHxA();
          tempXb = top->GetHxB();
          tempZa = top->GetHzA();
          tempZb = top->GetHxA();
      	}
      	else{
          tempXa.Add(&(top->GetHxA()));
          tempXb.Add(&(top->GetHxB()));
          tempZa.Add(&(top->GetHzA()));
          tempZb.Add(&(top->GetHzB()));
      	}
      }
      tempXa.SetName(Form("hXA%d",iGroup));
      tempXb.SetName(Form("hXB%d",iGroup));
      tempZa.SetName(Form("hZA%d",iGroup));
      tempZb.SetName(Form("hZB%d",iGroup));

    }
    fArrHisto.AddAt(new TH2F(tempXa),iGroup*4);
    fArrHisto.AddAt(new TH2F(tempXb),iGroup*4+1);
    fArrHisto.AddAt(new TH2F(tempZa),iGroup*4+2);
    fArrHisto.AddAt(new TH2F(tempZb),iGroup*4+3);
  }
  static TF1* gs = new TF1("gs","gaus",-50,50);
  static TF1* gs2 = new TF1("gs2","gaus",-50,50);
  printf("Taking data from histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    //X projection
    Int_t fitStatusX=0;
    TH2* hXA = (TH2F*)fArrHisto.At(iGroup*4);
    TH1* hdx = hXA->ProjectionY("hdx");
    gs->SetParameters(hdx->GetMaximum(),hdx->GetMean(),hdx->GetRMS());
    if((hdx->GetEntries())<100) fitStatusX = hdx->Fit("gs","ql");
    else fitStatusX = hdx->Fit("gs","q");
    //Z projection
    Int_t fitStatusZ=0;
    TH2* hZA = (TH2F*)fArrHisto.At(iGroup*4+2);
    TH1* hdz = hZA->ProjectionY("hdz");
    gs2->SetParameters(hdz->GetMaximum(),hdz->GetMean(),hdz->GetRMS());
    if((hdz->GetEntries())<100) fitStatusZ = hdz->Fit("gs2","ql");
    else fitStatusZ = hdz->Fit("gs2","q");
    //*******************************************
    for(Int_t i=0; i<nPatterns; i++){
      Topology* top = (Topology*)fArrTopologies.At(i);
      if(top->GetGroupID()==iGroup){
      	//
        top->SetHxA(*((TH2F*)fArrHisto.At(iGroup*4)));
        top->SetHxB(*((TH2F*)fArrHisto.At(iGroup*4+1)));
        top->SetHzA(*((TH2F*)fArrHisto.At(iGroup*4+2)));
        top->SetHzB(*((TH2F*)fArrHisto.At(iGroup*4+3)));
        //x chunk
      	if(fitStatusX==0){
      	  top->SetFitStuff(gs->GetParameter(1),Topology::kDeltaXmean);
      	  top->SetFitStuff(gs->GetParError(1),Topology::kDeltaXmeanErr);
      	  top->SetFitStuff(gs->GetParameter(2),Topology::kDeltaXsigma);
      	  top->SetFitStuff(gs->GetParError(2),Topology::kDeltaXsigmaErr);
      	  top->SetFitStuff(gs->GetChisquare(),Topology::kChi2x);
      	  Int_t varNDFx = gs->GetNDF();
      	  if(varNDFx>=0){
      	    top->SetFitStuff(varNDFx,Topology::kNDFx);
      	  }
      	  else{
      	    top->SetFitStuff(0.,Topology::kNDFx);
      	  }
      	}
      	else{
      	  top->SetFitStuff(0.,Topology::kDeltaXmean);
      	  top->SetFitStuff(0.,Topology::kDeltaXmeanErr);
      	  top->SetFitStuff(0.,Topology::kDeltaXsigma);
      	  top->SetFitStuff(0.,Topology::kDeltaXsigmaErr);
      	  top->SetFitStuff(-1,Topology::kChi2x);
      	}
      	//z chunk
      	if(fitStatusZ==0){
      	  top->SetFitStuff(gs2->GetParameter(1),Topology::kDeltaZmean);
      	  top->SetFitStuff(gs2->GetParError(1),Topology::kDeltaZmeanErr);
      	  top->SetFitStuff(gs2->GetParameter(2),Topology::kDeltaZsigma);
      	  top->SetFitStuff(gs2->GetParError(2),Topology::kDeltaZsigmaErr);
      	  top->SetFitStuff(gs2->GetChisquare(),Topology::kChi2z);
      	  Int_t varNDFz = gs2->GetNDF();
      	  if(varNDFz>=0){
      	    top->SetFitStuff(varNDFz,Topology::kNDFz);
      	  }
      	  else{
      	    top->SetFitStuff(0.,Topology::kNDFz);
      	  }
      	}
      	else{
      	  top->SetFitStuff(0.,Topology::kDeltaZmean);
      	  top->SetFitStuff(0.,Topology::kDeltaZmeanErr);
      	  top->SetFitStuff(0.,Topology::kDeltaZsigma);
      	  top->SetFitStuff(0.,Topology::kDeltaZsigmaErr);
      	  top->SetFitStuff(-1,Topology::kChi2z);
      	}
      }
    }
  }
}

Int_t TopDatabase::FromCluster2GroupID(const AliITSMFTClusterPix &cl) const{
  Topology top(cl);
  Int_t hashcode = top.GetHash();
  //Looking up in the database with intepolation search
  Int_t low = 0;
  Int_t high = fN-1;
  Int_t interIndex=-1;
  Int_t min=-2147483648;//minimum Int_t value
  Int_t max=2147483647;
  while(1){
    if(low>high || hashcode<min || hashcode>max){
      interIndex=-1;
      break;
    }
    Int_t guess;
    if(high==low) guess=high;
    else{
      Int_t size = high-low;
      Int_t offset = (Int_t)(((size-1)*((Long_t)hashcode-(Long_t)min))/((Long_t)max-(Long_t)min));
      guess = low+offset;
    }
    if( ( (Topology*)fArrTopologies.At(guess) )->GetHash()==hashcode){
      interIndex=guess;
      break;
    }
    else if( ( (Topology*)fArrTopologies.At(guess) )->GetHash()>hashcode){
      high = guess-1;
      if(high==-1) cout << "hashcode: " << hashcode << endl;
      max = ( (Topology*)fArrTopologies.At(guess-1) )->GetHash();
    }
    else{
      low = guess+1;
      min = ( (Topology*)fArrTopologies.At(guess+1) )->GetHash();
    }
  }
  if(interIndex==-1){//junk case
    return -1;
  }
  //Solving clashes if necessary
  if(( (Topology*)fArrTopologies.At(interIndex) )->GetFlag()==1){
    //printf("Clash found.\n");
    Bool_t IndexFound = kFALSE;
    Int_t guessIndex = interIndex;
    while(((Topology*)fArrTopologies.At(guessIndex))->GetHash()==hashcode && IndexFound==kFALSE){
      string str = ((Topology*)fArrTopologies.At(guessIndex))->GetPattern();
      if(top.GetPattern().compare(2,4,str,2,4)==0){
      	IndexFound = kTRUE;
      	interIndex = guessIndex;
      }
      guessIndex--;
    }
    guessIndex = interIndex;
    while(((Topology*)fArrTopologies.At(guessIndex))->GetHash()==hashcode && IndexFound==kFALSE){
      string str = ((Topology*)fArrTopologies.At(guessIndex))->GetPattern();
      if(top.GetPattern().compare(2,4,str,2,4)==0){
      	IndexFound = kTRUE;
      	interIndex = guessIndex;
      }
      guessIndex++;
    }
    if(IndexFound==kFALSE){
      printf("Index not found after clash\n");
      exit(1);
    }
  }
  return ( (Topology*)fArrTopologies.At(interIndex) )->GetGroupID();
}
