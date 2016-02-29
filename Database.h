#ifndef DATABASE_H
#define DATABASE_H
#include "TObject.h"
#include "TBits.h"
#include <Riostream.h>
#include "TArrayI.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "AliITSUClusterPix.h"

class Database :public TObject {

 public:
  Database();
  Database(const Database &ogg);
  ~Database();
  Database& operator=(const Database &ogg);
  void AccountTopology(const AliITSUClusterPix &cluster, Float_t dY, Float_t dZ, Float_t alpha, Float_t beta);  

  //________________getters___________________

  Int_t GetN() const {return fN;}
  TArrayI GetArrPattID() const {return fArrPattID;}
  TObjArray GetArrStore() const {return fArrStore;}
  TArrayF GetArrFreq() const {return fArrFreq;}
  TArrayF GetArrzCOGPix() const {return fArrzCOGPix;}
  TArrayF GetArrxCOGPix() const {return fArrxCOGPix;}
  TArrayF GetArrzCOGshift() const {return fArrzCOGshift;}
  TArrayF GetArrxCOGshift() const {return fArrxCOGshift;}
  TArrayI GetArrNpix() const {return fArrNpix;}
  TArrayI GetArrNcol() const {return fArrNcol;}
  TArrayI GetArrNrow() const {return fArrNrow;}
  TArrayI GetArrHash() const {return fArrHash;}
  TArrayF GetArrDeltaZmean() const {return fArrDeltaZmean;}
  TArrayF GetArrDeltaZmeanErr() const {return fArrDeltaZmeanErr;}
  TArrayF GetArrDeltaXmean() const {return fArrDeltaXmean;}
  TArrayF GetArrDeltaXmeanErr() const {return fArrDeltaXmeanErr;}
  TArrayF GetArrDeltaZsigma() const {return fArrDeltaZsigma;}
  TArrayF GetArrDeltaZsigmaErr() const {return fArrDeltaZsigmaErr;}
  TArrayF GetArrDeltaXsigma() const {return fArrDeltaXsigma;}
  TArrayF GetArrDeltaXsigmaErr() const {return fArrDeltaXsigmaErr;}
  TArrayF GetArrChi2x() const {return fArrChi2x;}
  TArrayF GetArrChi2z() const {return fArrChi2z;}
  TArrayI GetArrNDFx() const {return fArrNDFx;}
  TArrayI GetArrNDFz() const {return fArrNDFz;}
  TArrayI GetArrGroupID() const {return fArrGroupID;}
  TArrayI GetArrFlag() const {return fArrFlag;}
  TArrayI GetArrPartialTop() const {return fArrPartialTop;}
  TArrayI GetArrCount() const {return fArrCount;}
  Int_t GetTotClusters() const {return fTotClusters;}
  Float_t GetThreshold() const {return fThreshold;}
  Int_t GetOverThr() const {return fOverThr;}
  Int_t GetNGroups() const {return fNGroups;}
  Int_t GetNmax() const {return fNmax;}
  TObjArray GetArrHisto() const {return fArrHisto;}
  void LoadDB(const char* fname);
  void SaveDB(const char* fname);

  //*************************************************************

  void EndAndSort();//to end the database and sort key wrt hashes, in ascending order
  void PrintDB(const char* output = "Database.txt") const; //print the database on a txt file
  void SetThresholdCumulative(Float_t cumulative);
  //Threshold is the frequency for which you have a fraction = cumulative of topology not in groups
  void SetThreshold(Float_t thr);
  void Grouping(Int_t NumberofShiftXbins, Int_t NumberofShiftZbins);//return patterns over threshold
  void Grouping(Float_t threshold, Int_t NumberofShiftXbins=10, Int_t NumberofShiftZbins=10);
  void SetNmax(Int_t a) { fNmax = a;}

  //*************************************************************

  static void Top2Word(const TBits* top, UChar_t* Word);
  //Word: 1st byte = row span; 2nd = column span; others: pattern.
  //The length must be the minimum possible.
  static Int_t Top2Int(const TBits* top);//retunr the first 32 bits of the topology
  static void Word2Top(const UChar_t* Word, TBits &top);
  static std::ostream& printTop(TBits top, std::ostream &out);
  static Bool_t compTop(TBits top1, TBits top2);
  static UInt_t FuncMurmurHash2(const void * key, Int_t len);
  Int_t FromCluster2GroupID(const AliITSUClusterPix &cl) const;


 private:
  Int_t fN; //length of arrays
  TArrayI fArrPattID; //array ofpattern IDs according to frequency ( 0-> most frequent)
  TObjArray fArrStore; //array of patterns (TBits* format, unique ID = rs<<16 + cs)
  TArrayF fArrFreq; //array of frequencies
  TArrayF fArrzCOGPix; //z position of the centre of the pixel containing COG
  TArrayF fArrxCOGPix;
  TArrayF fArrzCOGshift; //z distance between COG and center of pixel containig it
  TArrayF fArrxCOGshift;
  TArrayI fArrNpix; //number of fired pixels
  TArrayI fArrNcol; //number of columns
  TArrayI fArrNrow;
  TArrayI fArrHash; //hashes corresponding to topologies
  TArrayF fArrDeltaZmean; //mean values of z distance between MC-truth and COG, from gaussian fit
  TArrayF fArrDeltaZmeanErr; //related error
  TArrayF fArrDeltaXmean;
  TArrayF fArrDeltaXmeanErr;
  TArrayF fArrDeltaZsigma; //sigma...
  TArrayF fArrDeltaZsigmaErr; //relative error...
  TArrayF fArrDeltaXsigma;
  TArrayF fArrDeltaXsigmaErr;
  TArrayF fArrChi2x; //Chi-square values of previous fits
  TArrayF fArrChi2z;
  TArrayI fArrNDFx; //Number of degrees of fredom of previous fits
  TArrayI fArrNDFz;
  TArrayI fArrGroupID; //array of group-IDs
  TArrayI fArrFlag; //0 if the hash is unique to that topology, 1 there is a clash
  TArrayI fArrPartialTop; //first 32 pixels of the topology
  Int_t fTotClusters; //number of clusters on wich DB has been constructed
  TArrayI fArrCount; //array of topologies counts
  Float_t fThreshold;//frequency threshold
  Int_t fOverThr;//number of patterns topologies over threshold
  Int_t fNGroups;//number of groups
  TObjArray fArrHisto;
  Int_t fNmax;//patterns above this number (included) belong to a "junk" bin
  void ExpandDB();

  ClassDef(Database,6)
};

#endif
