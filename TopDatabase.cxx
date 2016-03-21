#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include <limits>
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

using namespace std;

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

void TopDatabase::AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ, float alpha, float beta){
  Topology top(cluster);
  bool newPatt = true;
  bool Junk = false;
  fTotClusters++;
  int indTop = -1;
  for(int ip=0;ip<fN;ip++) {
    string pattOld = ((Topology*)fArrTopologies.At(ip))->GetPattern();
    if(top.GetPattern().length() == pattOld.length() && memcmp(top.GetPattern().data(),pattOld.data(),pattOld.length())==0){
      newPatt = false;
      //cout<< "old" << endl;
      indTop = ip;
      break;
    }
  }
  if(newPatt){
    //cout << "new" << endl;
    if(fN == fNmax){ //Junk bin
      Junk = true;
    }
    else {
      this->ExpandDB(top);
      indTop=fN-1;
    }
  }
  if(Junk==true) return;
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

void TopDatabase::SetThreshold(float thr){
  fThreshold=thr;
  int nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(int i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    float tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());
  int over=0;
  float fr=0.;
  for(int j=0; j<fN; j++){
    fr=arrFreq[sortIndex[j]];
    if(fr<thr) break;
    else over++;
  }
  fOverThr=over;
}

void TopDatabase::SetThresholdCumulative(float cumulative){
  if(cumulative<=0. || cumulative >=1.) cumulative = 0.99;
  float totFreq = 0.;
  int nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(int i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    float tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TArrayF provvFreq;
  provvFreq.Set(fN);
  TMath::Sort(fN,arrFreq.GetArray(),sortIndex.GetArray());
  for(int j=0; j<fN; j++){
    provvFreq[j]=arrFreq[sortIndex[j]];
  }
  int over=0;
  while(totFreq < cumulative){
    totFreq+=provvFreq[over++];
  }
  fThreshold=provvFreq[--over];
  while(provvFreq[over]==fThreshold) over--;
  fThreshold=provvFreq[over++];
  fOverThr = over;
}

void TopDatabase::EndAndSort(int mode){
  Topology::SetMode(mode);
  int nPatterns = fN;
  fArrTopologies.Sort();
  bool alreadyProcessed[nPatterns];
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(int i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    float tempFreq = ((float)(top->GetCounts()))/fTotClusters;
    top->SetFreq(tempFreq);
    arrFreq[i] = tempFreq;
    top->SetFlag(0);
    alreadyProcessed[i]=false;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());
  for(int i=0; i<nPatterns; i++){
    if(alreadyProcessed[i]==true) continue;
    else alreadyProcessed[i]=true;
    Topology* topFreq = (Topology*)fArrTopologies.At(sortIndex[i]);
    topFreq->SetPattID(i);
    Topology* top = (Topology*)fArrTopologies.At(i);
    unsigned long refHash = top->GetHash();
    for(int j=i+1; j<nPatterns; j++){
      Topology* top2comp = (Topology*)fArrTopologies.At(j);
      int counter = 0;
      if(top2comp->GetHash()==refHash){
      	alreadyProcessed[j]=true;
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
  for(int i=0; i<fN; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    o << "               POSITION " << i << endl;
    o << Form("\nPattID %4d R:%d C: %d Freq:%f Hash:%lu\n\n", top->GetPattID(),top->GetRowSpan(),top->GetColumnSpan(),top->GetFreq(),top->GetHash());
    top->printTop(o);
    o << "\nxCentre: " << top->GetxCOGPix() << " + " << top->GetxCOGshift()
      <<" zCentre: " << top->GetzCOGPix() << " + " << top->GetzCOGshift()
      <<" groupID: " << top->GetGroupID()<< "\n\n...............................................\n";
  }
  o.close();
}

void TopDatabase::Grouping(int NumberofShiftXbins, int NumberofShiftZbins){
  printf("\nGROUPING\n\n");
  float threshold=fThreshold;
  int nPatterns=fN;
  TArrayI shiftXID;
  TArrayI shiftZID;
  int notINgorups=0;
  float shiftXwidth=1./NumberofShiftXbins; //(fraction of the pitch)
  float shiftZwidth=1./NumberofShiftZbins; //(fraction of the pitch)
  float xOffset = shiftXwidth/2;
  float zOffset = shiftZwidth/2;

  shiftXID.Set(nPatterns);
  shiftZID.Set(nPatterns);

  TH1F xShiftDistro("xShiftDistro", "x position", NumberofShiftXbins+1 , -0.5-xOffset, 0.5+xOffset);
  xShiftDistro.SetDirectory(0);

  TH1F zShiftDistro("zShiftDistro", "z position", NumberofShiftZbins+1, -0.5-zOffset, 0.5+zOffset);
  zShiftDistro.SetDirectory(0);

  printf("Assigning shift-IDs:\n");
  for(int i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    xShiftDistro.Fill(top->GetxCOGshift());
    shiftXID[i]=xShiftDistro.FindBin(top->GetxCOGshift());
    zShiftDistro.Fill(top->GetzCOGshift());
    shiftZID[i]=zShiftDistro.FindBin(top->GetzCOGshift());
  }
  int tempgroupID=0;
  printf("Processing patterns over threshold:\n");
  //first assigning groupID to patterns whose frequency is above the treshold
  for(int t=0; t<nPatterns;t++){
    if(t%100==0)printf("%d / %d\n ", t, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(t);
    if(top->GetFreq()>=threshold){
      top->SetGroupID(tempgroupID);
      tempgroupID++;
      notINgorups++;
    }
    else break;
  }
  printf("Processing patterns under threshold:\n");
  for(int i=0; i<nPatterns; i++){
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(top->GetGroupID()!=-1) continue;
    top->SetGroupID(tempgroupID);
    //if(i%10000==0)printf("Assigning group ID %d / %d\n ", i, nPatterns);
    for(int j=i+1; j<nPatterns; j++){
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
  for(int iGroup=0; iGroup<fNGroups; iGroup++){
    TH2F tempXa(Form("hXA%d",iGroup),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempXb(Form("hXB%d",iGroup),"#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempZa(Form("hZA%d",iGroup),"#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
    TH2F tempZb(Form("hZB%d",iGroup),"#DeltaZ vs #beta",10,0,TMath::Pi()/2,50,-30,30);
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    bool FirstMatch = true;
    for(int i=0; i<nPatterns; i++){
      // printf("iGroup: %d i: %d\n", iGroup, i);
      Topology* top = (Topology*)fArrTopologies.At(i);
      if(top->GetGroupID()==iGroup){
      	if(FirstMatch==true){
      	  FirstMatch==false;
          tempXa = top->GetHxA();
          tempXb = top->GetHxB();
          tempZa = top->GetHzA();
          tempZb = top->GetHzB();
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
  for(int iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    //X projection
    int fitStatusX=0;
    TH2* hXA = (TH2F*)fArrHisto.At(iGroup*4);
    TH1* hdx = hXA->ProjectionY("hdx");
    gs->SetParameters(hdx->GetMaximum(),hdx->GetMean(),hdx->GetRMS());
    if((hdx->GetEntries())<100) fitStatusX = hdx->Fit("gs","ql");
    else fitStatusX = hdx->Fit("gs","q");
    //Z projection
    int fitStatusZ=0;
    TH2* hZA = (TH2F*)fArrHisto.At(iGroup*4+2);
    TH1* hdz = hZA->ProjectionY("hdz");
    gs2->SetParameters(hdz->GetMaximum(),hdz->GetMean(),hdz->GetRMS());
    if((hdz->GetEntries())<100) fitStatusZ = hdz->Fit("gs2","ql");
    else fitStatusZ = hdz->Fit("gs2","q");
    //*******************************************
    for(int i=0; i<nPatterns; i++){
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
      	  int varNDFx = gs->GetNDF();
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
      	  int varNDFz = gs2->GetNDF();
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

int TopDatabase::FromCluster2GroupID(const AliITSMFTClusterPix &cl) const{
  Topology top(cl);
  unsigned long hashcode = top.GetHash();
  //Looking up in the database with intepolation search
  int low = 0;
  int high = fN-1;
  int interIndex=-1;
  unsigned long min=0;//minimum unsigned long value
  unsigned long max=std::numeric_limits<unsigned long>::max();
  while(1){
    if(low>high || hashcode<min || hashcode>max){
      interIndex=-1;
      break;
    }
    int guess;
    if(high==low) guess=high;
    else{
      int size = high-low;
      int offset = (int)(((size-1)*((unsigned long)hashcode-(unsigned long)min))
      /((unsigned long)max-(unsigned long)min));
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
  string str = top.GetPattern();
  string str1 = ((Topology*)fArrTopologies.At(interIndex))->GetPattern();
  //if(str.length()!=str1.length() || str.compare(0,str1.length()-1,str1,0,str1.length()-1)!=0) cout << "Error during research" << endl;
  return ( (Topology*)fArrTopologies.At(interIndex) )->GetGroupID();
}

void TopDatabase::BuildMap(){
  fMap.clear();
  for(int i=0; i<fN; i++){
    string &str = ((Topology*)fArrTopologies.At(i))->GetPattern();
    unsigned long hash = (unsigned long)((Topology*)fArrTopologies.At(i))->GetHash();
    int groupID = ((Topology*)fArrTopologies.At(i))->GetGroupID();
    if(!fMap.insert(pair<unsigned long, int>(hash,groupID)).second) cout << "errore" << endl;
  }
}

int TopDatabase::FromCluster2GroupIDMap(const AliITSMFTClusterPix &cl) const{
  Topology top(cl);
  string &str = top.GetPattern();
  unsigned long hash = (unsigned long)top.GetHash();
  return fMap.find(hash)->second;
}

std::ostream& TopDatabase::showMap(std::ostream &out){
  typedef map<unsigned long, int>::const_iterator MapIterator;
  for(MapIterator iter = fMap.begin(); iter != fMap.end(); iter++){
    out << "Hash: " << iter->first << " GroupID: " << iter->second << endl;
  }
}

void TopDatabase::CompareMap(){
  int i=0;
  for(auto &p : fMap){
    unsigned long MapHash = p.first;
    unsigned long TopHash = ((Topology*)fArrTopologies.At(i))->GetHash();
    if(MapHash!=TopHash) printf("%d)  Different Values: %lu %lu\n",i,MapHash,TopHash);
    //cout << i << ")  Different Values: " << MapHash << " " << TopHash << endl;
    //else cout << i << ")  " << MapHash << " " << TopHash << endl;
    else printf("%d)  %lu %lu\n",i,MapHash,TopHash);
    ((Topology*)fArrTopologies.At(i))->printTop(cout);
    i++;
  }
}
