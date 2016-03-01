#ifndef TOPDATABASE_H
#define TOPDATABASE_H
#include "TObject.h"
#include "TBits.h"
#include "./Topology.h"
#include <Riostream.h>
#include "TArrayI.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "AliITSMFTClusterPix.h"
#include <map>

using namespace std;

class TopDatabase : public TObject {

 public:

  enum SortMode_t{kFrequency=0, kHashes=1};//fMode
  enum FitIndex_t{kDeltaZmean=0, kDeltaZmeanErr=1, kDeltaXmean=2, kDeltaXmeanErr=3, kDeltaZsigma=4, kDeltaZsigmaErr=5,
		  kDeltaXsigma=6, kDeltaXsigmaErr=7, kChi2z=8, kChi2x=9, kNDFx=10, kNDFz=11, kFitLength=12}; //position in Topology fArrFit

  TopDatabase();
  TopDatabase(TopDatabase &ogg);
  ~TopDatabase();
  void AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ, float alpha, float beta);

  int GetN() const {return fN;}
  int GetTotClusters() const {return fTotClusters;}
  float GetThreshold() const {return fThreshold;}
  int GetOverThr() const {return fOverThr;}
  int GetNGroups() const {return fNGroups;}
  int GetNmax() const {return fNmax;}

  void SetNmax(int a) { fNmax = a;}
  void SetThresholdCumulative(float cumulative);//Threshold is the frequency for which you have a fraction = cumulative of topology not in groups
  void SetThreshold(float thr);
  void EndAndSort(int mode = kHashes);//to end the database and sort key wrt hashes, in ascending order
  void PrintDB(const char* output = "Database.txt") const; //print the database on a txt file
  void Grouping(int NumberofShiftXbins, int NumberofShiftZbins);//return patterns over threshold
  void BuildMap();
  int FromCluster2GroupID(const AliITSMFTClusterPix &cl) const;
  int FromCluster2GroupIDMap(const AliITSMFTClusterPix &cl) const;
  std::ostream& showMap(std::ostream &out);
  void CompareMap();

private:
  string fPattern;
  int fN; //length of arrays
  TObjArray fArrTopologies;//array of topologies (class Topology)
  int fTotClusters;
  float fThreshold;//frequency threshold
  int fOverThr;//number of patterns topologies over threshold
  int fNGroups;
  int fNmax;//patterns above this number (included) belong to a "junk" bin
  map<unsigned long,int> fMap;
  TObjArray fArrHisto;
  TObjArray* GetArrTopologies() {return &fArrTopologies;}
  TObjArray* GetArrHisto() {return &fArrHisto;}
  void ExpandDB(const Topology &top);

ClassDef(TopDatabase,2)
};

#endif
