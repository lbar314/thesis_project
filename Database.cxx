#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Database.h"
#include "TString.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "AliITSUClusterPix.h"

ClassImp(Database)

Database::Database():TObject(),fN(0), fArrPattID(), fArrFreq(), fArrzCOGPix(), fArrxCOGPix(), fArrzCOGshift(), fArrxCOGshift(), fArrNpix(), fArrNcol(), fArrNrow(), fArrHash(), fArrStore(), fArrDeltaZmean(), fArrDeltaZmeanErr(), fArrDeltaXmean(), fArrDeltaXmeanErr(), fArrDeltaZsigma(), fArrDeltaZsigmaErr(), fArrDeltaXsigma(), fArrDeltaXsigmaErr(), fArrChi2x(), fArrChi2z(), fArrNDFx(), fArrNDFz(), fArrGroupID(), fArrFlag(), fArrPartialTop(), fTotClusters(0), fArrCount(), fThreshold(0.), fOverThr(0), fNGroups(0), fNmax(1e5) ,fArrHisto(){
}

Database::Database(const Database &ogg):TObject(){
  fN=ogg.GetN();
  fArrPattID=ogg.GetArrPattID();
  fArrFreq=ogg.GetArrFreq();
  fArrzCOGPix=ogg.GetArrzCOGPix();
  fArrxCOGPix=ogg.GetArrxCOGPix();
  fArrzCOGshift=ogg.GetArrzCOGshift();
  fArrxCOGshift=ogg.GetArrxCOGshift();
  fArrNpix=ogg.GetArrNpix();
  fArrNcol=ogg.GetArrNcol();
  fArrNrow=ogg.GetArrNrow();
  fArrHash=ogg.GetArrHash();
  fArrStore=ogg.GetArrStore();
  fArrDeltaZmean=ogg.GetArrDeltaZmean();
  fArrDeltaZmeanErr=ogg.GetArrDeltaZmeanErr();
  fArrDeltaXmean=ogg.GetArrDeltaXmean();
  fArrDeltaXmeanErr=ogg.GetArrDeltaXmeanErr();
  fArrDeltaZsigma=ogg.GetArrDeltaZsigma();
  fArrDeltaZsigmaErr=ogg.GetArrDeltaZsigmaErr();
  fArrDeltaXsigma=ogg.GetArrDeltaXsigma();
  fArrDeltaXsigmaErr=ogg.GetArrDeltaXsigmaErr();
  fArrChi2x=ogg.GetArrChi2x();
  fArrChi2z=ogg.GetArrChi2z();
  fArrNDFx=ogg.GetArrNDFx();
  fArrNDFz=ogg.GetArrNDFz();
  fArrGroupID=ogg.GetArrGroupID();
  fArrFlag=ogg.GetArrFlag();
  fArrPartialTop=ogg.GetArrPartialTop();
  fTotClusters=ogg.GetTotClusters();
  fArrCount=ogg.GetArrCount();
  fThreshold=ogg.GetThreshold();
  fOverThr=ogg.GetOverThr();
  fNGroups=ogg.GetNGroups();
  fArrHisto=ogg.GetArrHisto();
  fNmax=ogg.GetNmax();
}

void Database::ExpandDB(){
  fN++;
  fArrPattID.Set(fN);
  fArrStore.Expand(fN);
  fArrFreq.Set(fN);
  fArrzCOGPix.Set(fN);
  fArrxCOGPix.Set(fN);
  fArrzCOGshift.Set(fN);
  fArrxCOGshift.Set(fN);
  fArrNpix.Set(fN);
  fArrNcol.Set(fN);
  fArrNrow.Set(fN);
  fArrHash.Set(fN);
  fArrDeltaZmean.Set(fN);
  fArrDeltaZmeanErr.Set(fN);
  fArrDeltaXmean.Set(fN);
  fArrDeltaXmeanErr.Set(fN);
  fArrDeltaZsigma.Set(fN);
  fArrDeltaZsigmaErr.Set(fN);
  fArrDeltaXsigma.Set(fN);
  fArrDeltaXsigmaErr.Set(fN);
  fArrChi2x.Set(fN);
  fArrChi2z.Set(fN);
  fArrNDFx.Set(fN);
  fArrNDFz.Set(fN);
  fArrGroupID.Set(fN);
  fArrFlag.Set(fN);
  fArrPartialTop.Set(fN);
  fArrCount.Set(fN);
  fArrHisto.Expand(fN*4);

  //Creating histograms
  Int_t ind = fN-1;

  TH2* h0 = new TH2F(Form("hXalpha%d", ind),
		     Form("#DeltaX vs #alpha"),10,-1.1*TMath::Pi()/2, 1.1*TMath::Pi()/2,50,-30,30);
  h0->SetDirectory(0);
  h0->GetXaxis()->SetTitle("#alpha");
  h0->GetYaxis()->SetTitle("#DeltaX (#mum)");
  fArrHisto.AddAt(h0, ind*4);//This histogram is the first of the 4-histos series

  TH2* h1 = new TH2F(Form("hZalpha%d", ind),
		     Form("#DeltaZ vs #alpha"),10, -1.1* TMath::Pi()/2, 1.1* TMath::Pi()/2,50,-30,30);
  h1->SetDirectory(0);
  h1->GetXaxis()->SetTitle("#alpha");
  h1->GetYaxis()->SetTitle("#DeltaZ (#mum)");
  fArrHisto.AddAt(h1, ind*4+1);//This histogram is the second of the 4-histos series and so on...

  TH2* h2 = new TH2F(Form("hXbeta%d", ind),
		     Form("#DeltaX vs #beta"),10,-0.1, 1.1* TMath::Pi()/2,50,-30,30);
  h2->SetDirectory(0);
  h2->GetXaxis()->SetTitle("#beta");
  h2->GetYaxis()->SetTitle("#DeltaX (#mum)");
  fArrHisto.AddAt(h2, ind*4+2);

  TH2* h3 = new TH2F(Form("hZbeta%d", ind),
		     Form("#DeltaZ vs #beta"),10,-0.1, 1.1* TMath::Pi()/2,50,-30,30);
  h3->SetDirectory(0);
  h3->GetXaxis()->SetTitle("#beta");
  h3->GetYaxis()->SetTitle("#DeltaZ (#mum)");
  fArrHisto.AddAt(h3, ind*4+3);
}

Database::~Database(){
}

Database& Database::operator=(const Database &ogg){
  fN=ogg.GetN();
  fArrPattID=ogg.GetArrPattID();
  fArrFreq=ogg.GetArrFreq();
  fArrzCOGPix=ogg.GetArrzCOGPix();
  fArrxCOGPix=ogg.GetArrxCOGPix();
  fArrzCOGshift=ogg.GetArrzCOGshift();
  fArrxCOGshift=ogg.GetArrxCOGshift();
  fArrNpix=ogg.GetArrNpix();
  fArrNcol=ogg.GetArrNcol();
  fArrNrow=ogg.GetArrNrow();
  fArrHash=ogg.GetArrHash();
  fArrStore=ogg.GetArrStore();
  fArrDeltaZmean=ogg.GetArrDeltaZmean();
  fArrDeltaZmeanErr=ogg.GetArrDeltaZmeanErr();
  fArrDeltaXmean=ogg.GetArrDeltaXmean();
  fArrDeltaXmeanErr=ogg.GetArrDeltaXmeanErr();
  fArrDeltaZsigma=ogg.GetArrDeltaZsigma();
  fArrDeltaZsigmaErr=ogg.GetArrDeltaZsigmaErr();
  fArrDeltaXsigma=ogg.GetArrDeltaXsigma();
  fArrDeltaXsigmaErr=ogg.GetArrDeltaXsigmaErr();
  fArrChi2x=ogg.GetArrChi2x();
  fArrChi2z=ogg.GetArrChi2z();
  fArrNDFx=ogg.GetArrNDFx();
  fArrNDFz=ogg.GetArrNDFz();
  fArrGroupID=ogg.GetArrGroupID();
  fArrFlag=ogg.GetArrFlag();
  fArrPartialTop=ogg.GetArrPartialTop();
  fTotClusters=ogg.GetTotClusters();
  fArrCount=ogg.GetArrCount();
  fThreshold=ogg.GetThreshold();
  fOverThr=ogg.GetOverThr();
  fNGroups=ogg.GetNGroups();
  fArrHisto=ogg.GetArrHisto();
  fNmax=ogg.GetNmax();

  return *this;
}

void Database::AccountTopology(const AliITSUClusterPix &cluster, Float_t dY, Float_t dZ, Float_t alpha, Float_t beta){

  TBits Top;
  Top.Clear();
  Int_t rs = cluster.GetPatternRowSpan();
  Int_t cs = cluster.GetPatternColSpan();
  for(Int_t ir=0;ir<rs;ir++)
    for(Int_t ic=0;ic<cs;ic++)
      if(cluster.TestPixel(ir,ic)) Top.SetBitNumber(ir*cs+ic);
  Top.SetUniqueID((rs<<16)+cs);
  Bool_t newPatt = kTRUE;
  Bool_t Junk = kFALSE;
  fTotClusters++;
  Int_t indTop = -1;
  for(Int_t ip=0;ip<fN;ip++) {
    TBits pattOld = *(TBits*)fArrStore.At(ip);
    if(pattOld==Top && pattOld.GetUniqueID()==Top.GetUniqueID()){
      newPatt = kFALSE;
      indTop = ip;
      fArrCount[indTop]++;
      break;
    }
  }
  if(newPatt){
    if(fN == fNmax){ //Junk bin
      Junk = kTRUE;
    }
    else {
      TBits* pt = new TBits(Top);
      this->ExpandDB();
      indTop=fN-1;
      fArrStore.AddAt(pt,indTop);
      fArrCount[indTop]=1;
      //Defining Hash
      Int_t nBytes = rs*cs/8;
      if((rs*cs)%8 != 0) nBytes++;
      nBytes+=2; //first byte: number of rows; second byte: number of columns;
      const Int_t length = nBytes;
      UChar_t Word[length];
      for(int k=length; k--;) Word[k]=0;
      Top2Word(pt,Word);
      Int_t Hash = (Int_t)FuncMurmurHash2(Word,length);
      fArrHash[indTop]=Hash;
    }
  }
  if(Junk==kTRUE) return;
  TH2*  h0a = (TH2*)fArrHisto.At(indTop*4);
  h0a->Fill(alpha, dY);
  TH2*  h1a = (TH2*)fArrHisto.At(indTop*4+1);
  h1a->Fill(alpha, dZ);
  TH2*  h2a = (TH2*)fArrHisto.At(indTop*4+2);
  h2a->Fill(beta, dY);
  TH2*  h3a = (TH2*)fArrHisto.At(indTop*4+3);
  h3a->Fill(beta, dZ);
}

void Database::Top2Word(const TBits* top, UChar_t* Word){
  UInt_t UID = top->GetUniqueID();
  Int_t rs = UID>>16;
  Int_t cs = UID&0xffff;
  Word[0]=rs;
  Word[1]=cs;
  UChar_t tempChar=0;
  Int_t index=2;
  Int_t BitCounter=7;
  for(Int_t ir=0; ir<rs; ir++){
    for(Int_t ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	Word[index]=tempChar;
	tempChar=0;
	BitCounter=7;
	index++;
      }
      if(top->TestBitNumber(ir*cs+ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  Word[index]=tempChar;
}

void Database::Word2Top(const UChar_t* Word, TBits &top){
  Int_t rs = Word[0];
  Int_t cs = Word[1];
  Int_t nBytes = cs*rs/8;
  UChar_t tempChar=0;
  Int_t indBit=0;
  Int_t s=0;
  if((cs*rs)%8!=0) nBytes++;
  for(int i=2; i<nBytes+2; i++){
    tempChar=Word[i];
    s=128; //0b10000000
    while(s>0){
      if((tempChar&s)!=0) top.SetBitNumber(indBit);
      s/=2;
      indBit++;
    }
  }
  top.SetUniqueID((rs<<16)+cs);
}

Int_t Database::Top2Int(const TBits* top){
  Int_t output=0;
  for(Int_t i=0; i<32; i++){
    if(top->TestBitNumber(i)) output+=(1<<(31-i));
  }
  return output;
}

UInt_t Database::FuncMurmurHash2(const void* key, Int_t len){
  // 'm' and 'r' are mixing constants generated offline.
  const UInt_t m =0x5bd1e995;
  const Int_t r = 24;
  // Initialize the hash
  UInt_t h = 0;
  // Mix 4 bytes at a time into the hash
  const UChar_t* data = (const UChar_t *)key;
  //Int_t recIndex=0;
  while(len >= 4){
    UInt_t k = *(UInt_t*)data;
    k *= m;
    k ^= k >> r;
    k *= m;
    h *= m;
    h ^= k;
    data += 4;
    len -= 4;
  }
  // Handle the last few bytes of the input array
  switch(len){
  case 3: h ^= data[2] << 16;
  case 2: h ^= data[1] << 8;
  case 1: h ^= data[0];
    h *= m;
  };
  // Do a few final mixes of the hash to ensure the last few
  // bytes are well-incorporated.
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  return h;
}

Bool_t Database::compTop(TBits top1, TBits top2){
  if (top1 == top2 && top1.GetUniqueID()==top2.GetUniqueID()) return kTRUE;
  else return kFALSE;
}

std::ostream& Database::printTop(TBits top, std::ostream &out){
  UInt_t UID = top.GetUniqueID();
  Int_t rs = UID>>16;
  Int_t cs = UID&0xffff;
  for (Int_t ir=0;ir<rs;ir++){
    out << "|";
    for (Int_t ic=0; ic<cs; ic++) {
      out << Form("%c",top.TestBitNumber(ir*cs+ic) ? '+':' ');
    }
    out << ("|\n");
  }
  out<< endl;
}

void Database::EndAndSort(){
  //printf("EndAndSort. Number of entries in fArrHisto: %d", fArrHisto.GetEntries());
  TArrayI sortIndex;
  Int_t nPatterns = fN;
  sortIndex.Set(nPatterns);
  TMath::Sort(nPatterns, fArrHash.GetArray(), sortIndex.GetArray(), kFALSE);

  TObjArray tempArrStore = fArrStore;
  fArrStore.Clear();
  TObjArray tempArrHisto = fArrHisto;
  fArrHisto.Clear();
  TArrayI tempArrCount;
  tempArrCount.Set(nPatterns);
  TArrayI tempArrHash;
  tempArrHash.Set(nPatterns);

  for(Int_t i=0; i<nPatterns; i++){
    Int_t ind =sortIndex[i];
    fArrStore.AddAt(tempArrStore[ind],i);
    tempArrCount[i]=fArrCount[ind];
    tempArrHash[i]=fArrHash[ind];
    fArrHisto[i*4]=tempArrHisto[ind*4];
    fArrHisto[i*4+1]=tempArrHisto[ind*4+1];
    fArrHisto[i*4+2]=tempArrHisto[ind*4+2];
    fArrHisto[i*4+3]=tempArrHisto[ind*4+3];
    //Defining GOG position
    UInt_t idd = fArrStore[i]->GetUniqueID();
    Int_t firedPixels=0;
    Int_t tempxCOG=0;
    Int_t tempzCOG=0;
    Int_t rs = idd>>16;//It's the rows span
    Int_t cs = idd&0xffff;//It's the columns span
    for(Int_t ir=0;ir<rs;ir++){
      for(Int_t ic=0;ic<cs;ic++){
	if(((TBits*)fArrStore[i])->TestBitNumber(ir*cs+ic)){
	  firedPixels++;
	  tempxCOG+=ir;
	  tempzCOG+=ic;
	}
      }
    }
    Float_t xsh=Float_t((tempxCOG%firedPixels))/firedPixels; //distance between COG end centre of the pixel containing COG
    Float_t zsh=Float_t((tempzCOG%firedPixels))/firedPixels;
    tempxCOG/=firedPixels;
    tempzCOG/=firedPixels;
    if(xsh>0.5){
      tempxCOG+=1;
      xsh-=1;
    }
    if(zsh>0.5){
      tempzCOG+=1;
      zsh-=1;
    }
    Float_t xc=(Float_t) tempxCOG+0.5;
    Float_t zc=(Float_t) tempzCOG+0.5;
    fArrxCOGPix[i]=xc;
    fArrxCOGshift[i]=xsh;
    fArrzCOGPix[i]=zc;
    fArrzCOGshift[i]=zsh;
    fArrNpix[i]=firedPixels;
    fArrNrow[i]=rs;
    fArrNcol[i]=cs;
    fArrFreq[i]=(Float_t)tempArrCount[i]/fTotClusters;
  }
  TArrayI FreqIndex;
  FreqIndex.Set(nPatterns);
  TMath::Sort(nPatterns,fArrFreq.GetArray(),FreqIndex.GetArray());
  for(Int_t j=0; j<nPatterns; j++){
    fArrPattID[FreqIndex[j]]=j;
  }
  fArrCount=tempArrCount;
  fArrHash=tempArrHash;
  //Solving clashes dtoring the first 4 bytes of the topology
  Bool_t alreadyProcessed[nPatterns];
rr  for(Int_t l=0; l<nPatterns; l++){
    fArrFlag[l]=0;
    fArrPartialTop[l]=0;
    alreadyProcessed[l]=kFALSE;
  }
  for(Int_t i=0; i<nPatterns; i++){
    if(alreadyProcessed[i]==kTRUE) continue;
    else alreadyProcessed[i]=kTRUE;
    Int_t refHash = fArrHash[i];
    for(Int_t j=i+1; j<nPatterns; j++){
      Int_t counter = 0;
      if(fArrHash[j]==refHash){
	alreadyProcessed[j]=kTRUE;
	if(counter==0){
	  counter++;
	  TBits* top2study = (TBits*)fArrStore.At(i);
	  fArrFlag[i]=1;
	  fArrPartialTop[i]=Top2Int(top2study);
	}
	TBits* top2study = (TBits*)fArrStore.At(j);
	fArrFlag[j]=1;
	fArrPartialTop[j]=Top2Int(top2study);
      }
    }
  }
}

void Database::PrintDB(const char* output) const{
  ofstream o(output);
  o << "Number of topologies: " << fN << endl << endl;
  for(Int_t i=0; i<fN; i++){
    o << "               POSITION " << i << endl;
    o << Form("\nPattID %4d R:%d C: %d Freq:%f Hash:%d\n\n", (Int_t)fArrPattID[i],fArrNrow[i],fArrNcol[i],fArrFreq[i],(Int_t)fArrHash[i]);
    printTop(*(TBits*)fArrStore[i],o);
    o << "\nxCentre: " << fArrxCOGPix[i]  << " + " << fArrxCOGshift[i]
      <<" zCentre: " << fArrzCOGPix[i] << " + " << fArrzCOGshift[i]
      <<" groupID: " << fArrGroupID[i]<< "\n\n...............................................\n";
  }
  o.close();
}

void Database::SetThresholdCumulative(Float_t cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  Float_t totFreq = 0.;
  TArrayI indexArr;
  indexArr.Set(fN);
  TArrayF provvFreq;
  provvFreq.Set(fN);
  TMath::Sort(fN,fArrFreq.GetArray(),indexArr.GetArray());
  for(Int_t j=0; j<fN; j++){
    provvFreq[j]=fArrFreq[indexArr[j]];
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

void Database::Grouping(Float_t threshold, Int_t NumberofShiftXbins, Int_t NumberofShiftZbins){
  printf("\nGROUPING\n\n");
  fThreshold=threshold;
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
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    xShiftDistro.Fill(fArrxCOGshift[i]);
    shiftXID[i]=xShiftDistro.FindBin(fArrxCOGshift[i]);
    zShiftDistro.Fill(fArrzCOGshift[i]);
    shiftZID[i]=zShiftDistro.FindBin(fArrzCOGshift[i]);
  }
  for(Int_t ID=0; ID<nPatterns; ID++){
    fArrGroupID[ID]=-1;
  }
  Int_t tempgroupID=0;
  printf("Processing patterns over threshold:\n");
  //first assigning groupID to patterns whose frequency is above the treshold
  for(Int_t t=0; t<nPatterns;t++){
    if(t%100==0)printf("%d / %d\n ", t, nPatterns);
    if(fArrFreq[t]>=threshold){
      fArrGroupID[t]=tempgroupID;
      tempgroupID++;
      notINgorups++;
    }
  }
  printf("Processing patterns under threshold:\n");
  for(Int_t i=0; i<nPatterns; i++){
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    if(fArrGroupID[i]!=-1) continue;
    fArrGroupID[i]=tempgroupID;
    //if(i%10000==0)printf("Assigning group ID %d / %d\n ", i, nPatterns);
    for(Int_t j=i+1; j<nPatterns; j++){
      if(fArrGroupID[j]!=-1) continue;
      else {
	if(shiftXID[j]==shiftXID[i] && shiftZID[j]==shiftZID[i] &&
	   fArrNrow[i]==fArrNrow[j] && fArrNcol[i]==fArrNcol[j]) fArrGroupID[j]=tempgroupID;
      }
    }
    tempgroupID++;
  }
  fNGroups = tempgroupID;
  fOverThr = notINgorups;
  //**********************************************************//
  //                                                          //
  //                 OPERATIONS WITH HISTOGRAMS               //
  //                                                          //
  //**********************************************************//
  TObjArray NewHistoArr=fArrHisto;
  fArrHisto.Clear();

  printf("Number of groups: %d\n", fNGroups);

  printf("Summing group-histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    Bool_t FirstMatch = kTRUE;
    for(Int_t i=0; i<nPatterns; i++){
      // printf("iGroup: %d i: %d\n", iGroup, i);
      if(fArrGroupID[i]==iGroup){
	if(FirstMatch==kTRUE){
	  FirstMatch==kFALSE;
	  fArrHisto[iGroup*4]=NewHistoArr[i*4];
	  fArrHisto[iGroup*4+1]=NewHistoArr[i*4+1];
	  fArrHisto[iGroup*4+2]=NewHistoArr[i*4+2];
	  fArrHisto[iGroup*4+3]=NewHistoArr[i*4+3];
	}
	else{
	  ((TH2F*)fArrHisto[iGroup*4])->Add((TH2F*)NewHistoArr[i*4]);
	  ((TH2F*)fArrHisto[iGroup*4+1])->Add((TH2F*)NewHistoArr[i*4+1]);
	  ((TH2F*)fArrHisto[iGroup*4+2])->Add((TH2F*)NewHistoArr[i*4+2]);
	  ((TH2F*)fArrHisto[iGroup*4+3])->Add((TH2F*)NewHistoArr[i*4+3]);
	}
      }
    }
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
    TH2* hZA = (TH2F*)fArrHisto.At(iGroup*4+1);
    TH1* hdz = hZA->ProjectionY("hdz");
    gs2->SetParameters(hdz->GetMaximum(),hdz->GetMean(),hdz->GetRMS());
    if((hdz->GetEntries())<100) fitStatusZ = hdz->Fit("gs2","ql");
    else fitStatusZ = hdz->Fit("gs2","q");
    //*******************************************
    for(Int_t i=0; i<nPatterns; i++){
      if(fArrGroupID[i]==iGroup){
	//x chunk
	if(fitStatusX==0){
	  fArrDeltaXmean[i]= gs->GetParameter(1);
	  fArrDeltaXmeanErr[i]= gs->GetParError(1);
	  fArrDeltaXsigma[i]= gs->GetParameter(2);
	  fArrDeltaXsigmaErr[i]= gs->GetParError(2);
	  fArrChi2x[i] = gs->GetChisquare();
	  Int_t varNDFx = gs->GetNDF();
	  if(varNDFx>=0)
	    fArrNDFx[i] = varNDFx;
	  else{
	    fArrNDFx[i]=-1;
	  }
	}
	else{
	  fArrDeltaXmean[i]=0;
	  fArrDeltaXmeanErr[i]=0;
	  fArrDeltaXsigma[i]=0;
	  fArrDeltaXsigmaErr[i]=0;
	  fArrChi2x[i] = -1;
	}
	//z chunk
	if(fitStatusZ==0){
	  fArrDeltaZmean[i]= gs2->GetParameter(1);
	  fArrDeltaZmeanErr[i]= gs2->GetParError(1);
	  fArrDeltaZsigma[i]= gs2->GetParameter(2);
	  fArrDeltaZsigmaErr[i]= gs->GetParError(2);
	  fArrChi2z[i] = gs2->GetChisquare();
	  Int_t varNDFz = gs2->GetNDF();
	  if(varNDFz>=0)
	    fArrNDFz[i] = varNDFz;
	  else{
	    fArrNDFz[i]=-1;
	  }
	}
	else{
	  fArrDeltaZmean[i]=0;
	  fArrDeltaZmeanErr[i]=0;
	  fArrDeltaZsigma[i]=0;
	  fArrDeltaZsigmaErr[i]=0;
	  fArrChi2z[i] = -1;
	}
      }
    }
  }
}

void Database::Grouping(Int_t NumberofShiftXbins, Int_t NumberofShiftZbins){
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
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    xShiftDistro.Fill(fArrxCOGshift[i]);
    shiftXID[i]=xShiftDistro.FindBin(fArrxCOGshift[i]);
    zShiftDistro.Fill(fArrzCOGshift[i]);
    shiftZID[i]=zShiftDistro.FindBin(fArrzCOGshift[i]);
  }
  for(Int_t ID=0; ID<nPatterns; ID++){
    fArrGroupID[ID]=-1;
  }
  Int_t tempgroupID=0;
  printf("Processing patterns over threshold:\n");
  //first assigning groupID to patterns whose frequency is above the treshold
  for(Int_t t=0; t<nPatterns;t++){
    if(t%100==0)printf("%d / %d\n ", t, nPatterns);
    if(fArrFreq[t]>=threshold){
      fArrGroupID[t]=tempgroupID;
      tempgroupID++;
      notINgorups++;
    }
  }
  printf("Processing patterns under threshold:\n");
  for(Int_t i=0; i<nPatterns; i++){
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    if(fArrGroupID[i]!=-1) continue;
    fArrGroupID[i]=tempgroupID;
    //if(i%10000==0)printf("Assigning group ID %d / %d\n ", i, nPatterns);
    for(Int_t j=i+1; j<nPatterns; j++){
      if(fArrGroupID[j]!=-1) continue;
      else {
	if(shiftXID[j]==shiftXID[i] && shiftZID[j]==shiftZID[i] &&
	   fArrNrow[i]==fArrNrow[j] && fArrNcol[i]==fArrNcol[j]) fArrGroupID[j]=tempgroupID;
      }
    }
    tempgroupID++;
  }
  fNGroups = tempgroupID;
  fOverThr = notINgorups;
  //**********************************************************//
  //                                                          //
  //                 OPERATIONS WITH HISTOGRAMS               //
  //                                                          //
  //**********************************************************//
  TObjArray NewHistoArr=fArrHisto;
  fArrHisto.Clear();
  printf("Number of groups: %d\n", fNGroups);

  printf("Summing group-histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    Bool_t FirstMatch = kTRUE;
    for(Int_t i=0; i<nPatterns; i++){
      // printf("iGroup: %d i: %d\n", iGroup, i);
      if(fArrGroupID[i]==iGroup){
	if(FirstMatch==kTRUE){
	  FirstMatch==kFALSE;
	  fArrHisto[iGroup*4]=NewHistoArr[i*4];
	  fArrHisto[iGroup*4+1]=NewHistoArr[i*4+1];
	  fArrHisto[iGroup*4+2]=NewHistoArr[i*4+2];
	  fArrHisto[iGroup*4+3]=NewHistoArr[i*4+3];
	}
	else{
	  ((TH2F*)fArrHisto[iGroup*4])->Add((TH2F*)NewHistoArr[i*4]);
	  ((TH2F*)fArrHisto[iGroup*4+1])->Add((TH2F*)NewHistoArr[i*4+1]);
	  ((TH2F*)fArrHisto[iGroup*4+2])->Add((TH2F*)NewHistoArr[i*4+2]);
	  ((TH2F*)fArrHisto[iGroup*4+3])->Add((TH2F*)NewHistoArr[i*4+3]);
	}
      }
    }
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
    TH2* hZA = (TH2F*)fArrHisto.At(iGroup*4+1);
    TH1* hdz = hZA->ProjectionY("hdz");
    gs2->SetParameters(hdz->GetMaximum(),hdz->GetMean(),hdz->GetRMS());
    if((hdz->GetEntries())<100) fitStatusZ = hdz->Fit("gs2","ql");
    else fitStatusZ = hdz->Fit("gs2","q");
    //*******************************************
    for(Int_t i=0; i<nPatterns; i++){
      if(fArrGroupID[i]==iGroup){
	//x chunk
	if(fitStatusX==0){
	  fArrDeltaXmean[i]= gs->GetParameter(1);
	  fArrDeltaXmeanErr[i]= gs->GetParError(1);
	  fArrDeltaXsigma[i]= gs->GetParameter(2);
	  fArrDeltaXsigmaErr[i]= gs->GetParError(2);
	  fArrChi2x[i] = gs->GetChisquare();
	  Int_t varNDFx = gs->GetNDF();
	  if(varNDFx>=0)
	    fArrNDFx[i] = varNDFx;
	  else{
	    fArrNDFx[i]=-1;
	  }
	}
	else{
	  fArrDeltaXmean[i]=0;
	  fArrDeltaXmeanErr[i]=0;
	  fArrDeltaXsigma[i]=0;
	  fArrDeltaXsigmaErr[i]=0;
	  fArrChi2x[i] = -1;
	}
	//z chunk
	if(fitStatusZ==0){
	  fArrDeltaZmean[i]= gs2->GetParameter(1);
	  fArrDeltaZmeanErr[i]= gs2->GetParError(1);
	  fArrDeltaZsigma[i]= gs2->GetParameter(2);
	  fArrDeltaZsigmaErr[i]= gs->GetParError(2);
	  fArrChi2z[i] = gs2->GetChisquare();
	  Int_t varNDFz = gs2->GetNDF();
	  if(varNDFz>=0)
	    fArrNDFz[i] = varNDFz;
	  else{
	    fArrNDFz[i]=-1;
	  }
	}
	else{
	  fArrDeltaZmean[i]=0;
	  fArrDeltaZmeanErr[i]=0;
	  fArrDeltaZsigma[i]=0;
	  fArrDeltaZsigmaErr[i]=0;
	  fArrChi2z[i] = -1;
	}
      }
    }
  }
}

Int_t Database::FromCluster2GroupID(const AliITSUClusterPix &cl) const{

  //Passing from cluster to UChar_t*
  Int_t rs = cl.GetPatternRowSpan();
  Int_t cs = cl.GetPatternColSpan();
  Int_t nBytes = rs*cs/8;
  if((rs*cs)%8 != 0) nBytes++;
  nBytes+=2; //first byte: number of rows; second byte: number of columns;
  const Int_t length = nBytes;
  UChar_t Word[length];
  for(int k=length; k--;) Word[k]=0;
  Word[0]=rs;
  Word[1]=cs;
  UChar_t tempChar=0;
  Int_t index=2;
  Int_t BitCounter=7;
  for(Int_t ir=0; ir<rs; ir++){
    for(Int_t ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	Word[index]=tempChar;
	tempChar=0;
	BitCounter=7;
	index++;
      }
      if(cl.TestPixel(ir,ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  Word[index]=tempChar;
  //Creating Hash
  Int_t hashcode = (Int_t)FuncMurmurHash2(Word,length);
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
    if(fArrHash[guess]==hashcode){
      interIndex=guess;
      break;
    }
    else if(fArrHash[guess]>hashcode){
      high = guess-1;
      max = fArrHash[guess-1];
    }
    else{
      low = guess+1;
      min = fArrHash[guess+1];
    }
  }
  if(interIndex==-1){//junk case
    return -1;
  }
  //Solving clashes if necessary
  if(fArrFlag[interIndex]==1){
    //printf("Clash found.\n");
    Int_t Part = 0;
    Int_t BitPosition = 0;
    for(Int_t ir=0; ir<rs;ir++){
      for(Int_t ic=0; ic<cs; ic++){
	if(cl.TestPixel(ir,ic)) Part+=(1<<(31-BitPosition));
	BitPosition++;
	if(BitPosition==32) break;
      }
      if(BitPosition==32) break;
    }
    Bool_t IndexFound = kFALSE;
    Int_t guessIndex = interIndex;
    while(fArrHash[guessIndex]==hashcode && IndexFound==kFALSE){
      if(Part==fArrPartialTop[guessIndex]){
	IndexFound = kTRUE;
	interIndex = guessIndex;
      }
      guessIndex--;
    }
    guessIndex = interIndex;
    while(fArrHash[guessIndex]==hashcode && IndexFound==kFALSE){
      if(Part==fArrPartialTop[guessIndex]){
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
  return fArrGroupID[interIndex];
}

void Database::SetThreshold(Float_t thr){
  fThreshold=thr;
  TArrayI indexArr;
  indexArr.Set(fN);
  TArrayF provvFreq;
  provvFreq.Set(fN);
  TMath::Sort(fN,fArrFreq.GetArray(),indexArr.GetArray());
  Int_t over=0;
  for(Int_t j=0; j<fN; j++){
    provvFreq[j]=fArrFreq[indexArr[j]];
    if(provvFreq[j]<thr) break;
    else over++;
  }
  fOverThr=over;
}

void Database::LoadDB(const char* fname){
  TFile* fl = TFile::Open(fname);
  *this = *(Database*)fl->Get("DB");
  fl->Close();
}

void Database::SaveDB(const char* fname){
  TFile* flDB = TFile::Open(fname, "recreate");
  flDB->WriteObject(this,"DB","kSingleKey");
  flDB->Close();
}
