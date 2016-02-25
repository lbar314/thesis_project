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

int Topology::fMode = Topology::kHashes;

Topology::Topology():TObject(), fPattern(), fFiredPixels(0),fxCOGPix(0.), fzCOGPix(0.), fxCOGshift(0.), fzCOGshift(0.), fHash(0), fFreq(0.), fCounts(0),fGroupID(-1), fHxA(), fHxB(), fHzA(), fHzB(), fFlag(0), fPattID(-1){
  for(int i=0; i<kFitLength; i++) fArrFit[i]=0;
}

Topology::~Topology(){
}


Topology::Topology(const AliITSMFTClusterPix &cluster, int ID):TObject()
, fHxA(Form("hXA%d",ID),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30)
, fHxB(Form("hXB%d",ID),"#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30)
, fHzA(Form("hZA%d",ID),"#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30)
, fHzB(Form("hZB%d",ID),"#DeltaZ vs beta",10,0,TMath::Pi()/2,50,-30,30){
  int rs = cluster.GetPatternRowSpan();
  int cs = cluster.GetPatternColSpan();
	int tempxCOG = 0;
  int tempzCOG = 0;
  int tempFiredPixels = 0;
	//______________________________creating fPattern
	fPattern.push_back(rs);
	fPattern.push_back(cs);
  UChar_t tempChar=0;
  int BitCounter=7;
  for(int ir=0; ir<rs; ir++){
    for(int ic=0; ic<cs; ic++){
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
  float xsh=float((tempxCOG%tempFiredPixels))/tempFiredPixels; //distance between COG end centre of the pixel containing COG
  float zsh=float((tempzCOG%tempFiredPixels))/tempFiredPixels;
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
  fxCOGPix = (float) tempxCOG+0.5;
  fzCOGPix = (float) tempzCOG+0.5;
  fxCOGshift = xsh;
  fzCOGshift = zsh;
  fFiredPixels = tempFiredPixels;
  //__________________________________________________________Creating hash
  fHash=FuncMurmurHash2(fPattern.data(),fPattern.length());
  fFreq=0.; //WARNING: it is to set in a second time
  fCounts=0; //WARNING: it is to set in a second time
  fGroupID=-1; //WARNING: it is to set in a second time
	//_______________________________________________________Setting histograms
  fHxA.SetDirectory(0);
  fHxA.GetXaxis()->SetTitle("#alpha");
  fHxA.GetYaxis()->SetTitle("#DeltaX (#mum)");
  fHxB.SetDirectory(0);
  fHxB.GetXaxis()->SetTitle("#alpha");
  fHxB.GetYaxis()->SetTitle("#DeltaX (#mum)");
  fHzA.SetDirectory(0);
  fHzA.GetXaxis()->SetTitle("#alpha");
  fHzA.GetYaxis()->SetTitle("#DeltaZ (#mum)");
  fHzB.SetDirectory(0);
  fHzB.GetXaxis()->SetTitle("#beta");
  fHzB.GetYaxis()->SetTitle("#DeltaZ (#mum)");
  for(int i=0; i<kFitLength; i++) fArrFit[i]=0;
  fFlag=0;
  fPattID = -1;
}

Topology::Topology(const Topology &topo, int ID):TObject()
, fHxA(Form("hXA%d", ID),"#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30)
, fHxB(Form("hXB%d",ID),"#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30)
, fHzA(Form("hZA%d",ID),"#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30)
, fHzB(Form("hZB%d",ID),"#DeltaZ vs beta",10,0,TMath::Pi()/2,50,-30,30){
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
  for(int i=0; i<kFitLength; i++) fArrFit[i]=topo.GetFitStuff(i);
  fFlag=0;
  fPattID=-1;
}

std::ostream& Topology::printTop(std::ostream &out){
	int rs = fPattern[0];
	int cs = fPattern[1];
	out << "rs: " << rs << " cs: " << cs << " control: " << fPattern.length() << endl;
	UChar_t tempChar = 0;
	int s=0;
	int ic = 0;
  for (int i=2; i<fPattern.length(); i++){
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
	int rs = cluster.GetPatternRowSpan();
	int cs = cluster.GetPatternColSpan();
	for (int ir=0;ir<rs;ir++){
		out << "|";
		for (int ic=0; ic<cs; ic++) {
			out << Form("%c",cluster.TestPixel(ir,ic) ? '+':' ');
		}
		out << ("|\n");
	}
	out<< endl;
}

unsigned int Topology::FuncMurmurHash2(const void* key, int len){
  // 'm' and 'r' are mixing constants generated offline.
  const unsigned int m =0x5bd1e995;
  const int r = 24;
  // Initialize the hash
  unsigned int h = len^0xdeadbeef;
  // Mix 4 bytes at a time into the hash
  const UChar_t* data = (const UChar_t *)key;
  //int recIndex=0;
  while(len >= 4){
    unsigned int k = *(unsigned int*)data;
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


int Topology::Compare(const TObject* obj) const{
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

bool Topology::IsEqual(const TObject* obj) const{
  const Topology* top = (const Topology*)obj;
  if(fMode==kFrequency){
    if(fFreq == top->GetFreq()) return true;
    else return kFALSE;
  }
  if(fMode==kHashes){
    if(fHash == top->GetHash()) return true;
    else return kFALSE;
  }
  AliFatal(Form("Unknown mode for sorting: %d",fMode));
  return kFALSE;
}
