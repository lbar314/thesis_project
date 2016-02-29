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
  enum FitIndex_t{kDeltaZmean=0, kDeltaZmeanErr=1, kDeltaYmean=2, kDeltaYmeanErr=3, kDeltaZsigma=4, kDeltaZsigmaErr=5,
		  kDeltaYsigma=6, kDeltaYsigmaErr=7, kChi2z=8, kChi2y=9, kNDFy=10, kNDFz=11, kFitLength=12}; //position in fArrFit

  Topology();
  virtual ~Topology();
  Topology(const AliITSMFTClusterPix &cluster, int ID=-1);
  Topology(const Topology &topo, int ID=-1);

  bool IsEqual(const TObject* obj) const;
  bool IsSortable() const {return true;}
  int Compare(const TObject* obj) const;
  string& GetPattern() {return fPattern;};
  string  GetPattern() const {return fPattern;};
  int GetRowSpan() const {return fPattern[0];}
  int GetColumnSpan() const {return fPattern[1];}
  unsigned long GetHash() const {return fHash;}
  float GetFreq() const {return fFreq;}
  int GetCounts() const {return fCounts;}
  int GetGroupID() const {return fGroupID;}
  int GetFiredPixels() const {return fFiredPixels;}
  float GetyCOGPix() const {return fyCOGPix;}
  float GetzCOGPix() const {return fzCOGPix;}
  float GetyCOGshift() const {return fyCOGshift;}
  float GetzCOGshift() const {return fzCOGshift;}
  int GetMode() const {return fMode;}
  TH2F& GetHyA() {return fHyA;}
  TH2F& GetHyB() {return fHyB;}
  TH2F& GetHzA() {return fHzA;}
  TH2F& GetHzB() {return fHzB;}
  float GetFitStuff(int ind) const {return fArrFit[ind];}
  int GetFlag() const {return fFlag;}
  int GetPattID() const {return fPattID;}

  std::ostream& printTop(std::ostream &out);
  static unsigned int FuncMurmurHash2(const void * key, int len);
  static std::ostream& printCluster(const AliITSMFTClusterPix &cluster,std::ostream &out);

  void SetGroupID(int num){fGroupID=num;}
  void SetFreq(float num){fFreq=num;}
  void SetFitStuff(float value, int ind) {fArrFit[ind]=value;}
  void SetFlag(int num) {fFlag=num;}
  void IncreaseCounts(){fCounts++;}
  void SetHyA(TH2F &hist){fHyA=hist;}
  void SetHyB(TH2F &hist){fHyB=hist;}
  void SetHzA(TH2F &hist){fHzA=hist;}
  void SetHzB(TH2F &hist){fHzB=hist;}
  void SetPattID(int num) {fPattID=num;}
  static void SetMode(int mode) {fMode=mode;}

 private:

  string fPattern;
  int fFiredPixels;
  float fyCOGPix;
  float fzCOGPix;
  float fyCOGshift;
  float fzCOGshift;
  unsigned int fHash;
  int fCounts;
  float fFreq;
  int fGroupID;
  static int fMode; //DEFAULT kHashes
  TH2F fHyA;
  TH2F fHyB;
  TH2F fHzA;
  TH2F fHzB;
  float fArrFit[kFitLength];
  int fFlag;
  int fPattID;

  ClassDef(Topology,2)

};

#endif
