#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include "TObject.h"
#include "TBits.h"
#include "AliITSMFTClusterPix.h"
#include "TH2F.h"
#include "TH1F.h"
#include <Riostream.h>
#include <string>

class Topology :public TObject {

 public:
  enum SortMode_t{kFrequency=0, kHashes=1};//fMode
  enum FitIndex_t{kDeltaZmean=0, kDeltaZmeanErr=1, kDeltaXmean=2, kDeltaXmeanErr=3, kDeltaZsigma=4, kDeltaZsigmaErr=5,
		  kDeltaXsigma=6, kDeltaXsigmaErr=7, kChi2z=8, kChi2x=9, kNDFx=10, kNDFz=11, kFitLength=12}; //position in fArrFit

  Topology();
  virtual ~Topology();
  Topology(const AliITSMFTClusterPix &cluster);
  Topology(const Topology &topo);

  Bool_t IsEqual(const TObject* obj) const;
  Bool_t IsSortable() const {return kTRUE;}
  Int_t Compare(const TObject* obj) const;
  string GetPattern() const {return fPattern;};
  Int_t GetRowSpan() const {return fPattern[0];}
  Int_t GetColumnSpan() const {return fPattern[1];}
  Int_t GetHash() const {return fHash;}
  Float_t GetFreq() const {return fFreq;}
  Int_t GetCounts() const {return fCounts;}
  Int_t GetGroupID() const {return fGroupID;}
  Int_t GetFiredPixels() const {return fFiredPixels;}
  Float_t GetxCOGPix() const {return fxCOGPix;}
  Float_t GetzCOGPix() const {return fzCOGPix;}
  Float_t GetxCOGshift() const {return fxCOGshift;}
  Float_t GetzCOGshift() const {return fzCOGshift;}
  Int_t GetMode() const {return fMode;}
  TH2F* GetHxA() const {return fHxA;}
  TH2F* GetHxB() const {return fHxB;}
  TH2F* GetHzA() const {return fHzA;}
  TH2F* GetHzB() const {return fHzB;}
  Float_t GetFitStuff(Int_t ind) const {return fArrFit[ind];}
  Int_t GetFlag() const {return fFlag;}
  Int_t GetPattID() const {return fPattID;}

  std::ostream& printTop(std::ostream &out);
  static UInt_t FuncMurmurHash2(const void * key, Int_t len);
  static std::ostream& printCluster(const AliITSMFTClusterPix &cluster,std::ostream &out);

  void SetGroupID(Int_t num){fGroupID=num;}
  void SetFreq(Float_t num){fFreq=num;}
  void SetFitStuff(Float_t value, Int_t ind) {fArrFit[ind]=value;}
  void SetFlag(Int_t num) {fFlag=num;}
  void IncreaseCounts(){fCounts++;}
  void SetPattID(Int_t num) {fPattID=num;}
  static void SetMode(Int_t mode) {fMode=mode;}
  void DeleteHistos();
  void SetHxA(TH2F* ptr) {fHxA=ptr;}
  void SetHzA(TH2F* ptr) {fHzA=ptr;}
  void SetHxB(TH2F* ptr) {fHxB=ptr;}
  void SetHzB(TH2F* ptr) {fHzB=ptr;}

 private:

  string fPattern;
  Int_t fFiredPixels;
  Float_t fxCOGPix;
  Float_t fzCOGPix;
  Float_t fxCOGshift;
  Float_t fzCOGshift;
  Int_t fHash;
  Int_t fCounts;
  Float_t fFreq;
  Int_t fGroupID;
  static Int_t fMode; //DEFAULT kHashes
  TH2F* fHxA;
  TH2F* fHxB;
  TH2F* fHzA;
  TH2F* fHzB;
  Float_t fArrFit[kFitLength];
  Int_t fFlag;
  Int_t fPattID;

  ClassDef(Topology,1)

};

#endif
