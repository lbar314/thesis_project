#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "TBits.h"
#include "AliITSUClusterPix.h"
#include "TH2F.h"
#include "TH1F.h"
#include "Topology.h"
#include "AliLog.h"
#include "TString.h"

ClassImp(Topology)

Int_t Topology::fMode = Topology::kHashes;

Topology::Topology():TObject(), fPattern(), fFiredPixels(0),fxCOGPix(0.), fzCOGPix(0.), fxCOGshift(0.), fzCOGshift(0.), fHash(0), fFreq(0.), fCounts(0),fGroupID(-1), fHxA(0), fHxB(0), fHzA(0), fHzB(0), fFlag(0), fPattID(-1){
  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=0;
}

Topology::~Topology(){
  delete fHxA;
  delete fHxB;
  delete fHzA;
  delete fHzB;
}


Topology::Topology(const AliITSMFTClusterPix &cluster):TObject(){
  Int_t rs = cluster.GetPatternRowSpan();
  Int_t cs = cluster.GetPatternColSpan();
	Int_t tempxCOG = 0;
  Int_t tempzCOG = 0;
  Int_t tempFiredPixels = 0;
	//______________________________creating fPattern
	fPattern.push_back(rs);
	fPattern.push_back(cs);
  UChar_t tempChar=0;
  Int_t BitCounter=7;
  for(Int_t ir=0; ir<rs; ir++){
    for(Int_t ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	      fPattern.push_back(tempChar);
	      tempChar=0;
	      BitCounter=7;
      }
      if(cluster.TestPixel(ir,ic)){
				tempChar+=(1<<BitCounter);
				tempFiredPixels++;
				tempxCOG+=ir;
				tempzCOG+=ic;
			}
      BitCounter--;
    }
  }
	fPattern.push_back(tempChar);
  Float_t xsh=Float_t((tempxCOG%tempFiredPixels))/tempFiredPixels; //distance between COG end centre of the pixel containing COG
  Float_t zsh=Float_t((tempzCOG%tempFiredPixels))/tempFiredPixels;
  tempxCOG/=tempFiredPixels;
  tempzCOG/=tempFiredPixels;
  if(xsh>0.5){
    tempxCOG+=1;
    xsh-=1;
  }
  if(zsh>0.5){
    tempzCOG+=1;
    zsh-=1;
  }
  fxCOGPix = (Float_t) tempxCOG+0.5;
  fzCOGPix = (Float_t) tempzCOG+0.5;
  fxCOGshift = xsh;
  fzCOGshift = zsh;
  fFiredPixels = tempFiredPixels;
  //__________________________________________________________Creating hash
  fHash=(Int_t)FuncMurmurHash2(fPattern.data(),fPattern.length());
  fFreq=0.; //WARNING: it is to set in a second time
  fCounts=0; //WARNING: it is to set in a second time
  fGroupID=-1; //WARNING: it is to set in a second time
	//_______________________________________________________Creating histograms
  fHxA = new TH2F("hXA","#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHxA->SetDirectory(0);
  fHxA->GetXaxis()->SetTitle("#alpha");
  fHxA->GetYaxis()->SetTitle("#DeltaX (#mum)");

  fHxB = new TH2F("hXB","#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30);
  fHxB->SetDirectory(0);
  fHxB->GetXaxis()->SetTitle("#alpha");
  fHxB->GetYaxis()->SetTitle("#DeltaX (#mum)");

  fHzA = new TH2F("hZA","#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHzA->SetDirectory(0);
  fHzA->GetXaxis()->SetTitle("#alpha");
  fHzA->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  fHzB = new TH2F("hZB","#DeltaZ vs beta",10,0,TMath::Pi()/2,50,-30,30);
  fHzB->SetDirectory(0);
  fHzB->GetXaxis()->SetTitle("#beta");
  fHzB->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=0;
  fFlag=0;
  fPattID = -1;
}

Topology::Topology(const Topology &topo):TObject(){
	fPattern = topo.GetPattern();
  fFreq = topo.GetFreq();
  fCounts = topo.GetCounts();
  fHash = topo.GetHash();
  fGroupID = topo.GetGroupID();
  fFiredPixels = topo.GetFiredPixels();
  fxCOGPix = topo.GetxCOGPix();
  fzCOGPix = topo.GetzCOGPix();
  fxCOGshift = topo.GetxCOGshift();
  fzCOGshift = topo.GetzCOGshift();
  fHxA = new TH2F(*topo.GetHxA());
  fHxB = new TH2F(*topo.GetHxB());
  fHzA = new TH2F(*topo.GetHzA());
  fHzB = new TH2F(*topo.GetHzB());
  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=topo.GetFitStuff(i);
  fFlag=0;
  fPattID=-1;
}

std::ostream& Topology::printTop(std::ostream &out){

	Int_t rs = fPattern[0];
	Int_t cs = fPattern[1];
	out << "rs: " << rs << " cs: " << cs << " control: " << fPattern.length() << endl;
	UChar_t tempChar = 0;
	Int_t s=0;
	Int_t ic = 0;
  for (Int_t i=2; i<fPattern.length(); i++){
		tempChar = fPattern[i];
		s=128; //0b10000000
    while(s>0){
			if(ic%cs==0) out << "|";
			ic++;
      out << Form("%c", (tempChar&s)!=0 ? '+':' ');
      s/=2;
			if(ic%cs==0) out << "|" << endl;
			if(ic==(rs*cs)) break;
    }
		if(ic==(rs*cs)) break;
  }
  out<< endl;
}

std::ostream& Topology::printCluster(const AliITSMFTClusterPix &cluster, std::ostream &out){
	Int_t rs = cluster.GetPatternRowSpan();
	Int_t cs = cluster.GetPatternColSpan();
	for (Int_t ir=0;ir<rs;ir++){
		out << "|";
		for (Int_t ic=0; ic<cs; ic++) {
			out << Form("%c",cluster.TestPixel(ir,ic) ? '+':' ');
		}
		out << ("|\n");
	}
	out<< endl;
}

UInt_t Topology::FuncMurmurHash2(const void* key, Int_t len){
  // 'm' and 'r' are mixing constants generated offline.
  const UInt_t m =0x5bd1e995;
  const Int_t r = 24;
  // Initialize the hash
  UInt_t h = len^0xdeadbeef;
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


Int_t Topology::Compare(const TObject* obj) const{
  //Since Sort method of TObjArray sorts object in ascending order ( if +1 means higher and -1 lower),
  //it is necessary to invert the frequency outputs
  const Topology* top = (const Topology*)obj;
  if(fMode==kFrequency){
    if(fFreq < top->GetFreq()) return +1;
    else if(fFreq == top->GetFreq()) return 0;
    else return -1;
  }
  if(fMode==kHashes){
    if(fHash < top->GetHash()) return -1;
    else if(fHash == top->GetHash()) return 0;
    else return +1;
  }
  AliFatal(Form("Unknown mode for sorting: %d",fMode));
  return 0;
}

Bool_t Topology::IsEqual(const TObject* obj) const{
  const Topology* top = (const Topology*)obj;
  if(fMode==kFrequency){
    if(fFreq == top->GetFreq()) return kTRUE;
    else return kFALSE;
  }
  if(fMode==kHashes){
    if(fHash == top->GetHash()) return kTRUE;
    else return kFALSE;
  }
  AliFatal(Form("Unknown mode for sorting: %d",fMode));
  return kFALSE;
}

void Topology::DeleteHistos(){
  delete fHxA;
  delete fHxB;
  delete fHzA;
  delete fHzB;
}
