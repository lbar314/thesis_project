#ifndef BUILDDICTIONARY_H
#define BUILDDICTIONARY_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TH1F.h"
#include <map>
#include "AliITSMFTClusterPix.h"
#include "./MinimTopology.h"
#include "./Dictionary.h"

#define _STUDY_
#define _HISTO_ //in order to have a histogram with the ditribution of groupIDs

#ifdef _STUDY_
  struct TopologyInfo{
    int sizeX;
    int sizeZ;
    float xCOG;
    float zCOG;
    TH1F hdX;
    TH1F hdZ;
    int nPixels;
  };
#endif

struct groupTmp{
  unsigned long GrCounts;
  int tempGroupID;
  vector<unsigned long> GrHashes;
  //float tempErrX;
  //float tempErrZ;
  int groupID;
};

class BuildDictionary {

  public:
    #ifdef _HISTO_
      TH1F fHdist; //Distribution of groupIDs
    #endif

    #ifndef _STUDY_
      void AccountTopology(const AliITSMFTClusterPix &cluster);
    #else
      void AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ);
    #endif

    void SetThresholdCumulative(double cumulative);
    void Grouping();
    std::ostream& showMap(std::ostream &out);
    void PrintDictionary(string fname);

    int GetTotClusters() const {return fTotClusters;}
    int GetNotInGroups() const {return fNotInGroups;}
    int GetNGroups() const {return fNGroups;}

  private:
    map<unsigned long,pair<MinimTopology,unsigned long>> fMapTop; //<hash,<topology,counts>>,
    vector <pair<unsigned long,unsigned long>> fTopFreq; //<freq,hash>, needed to define threshold
    int fTotClusters;
    int fNGroups;
    int fNotInGroups;
    double fThreshold;
    Dictionary fDict;
    #ifdef _STUDY_
      map<long unsigned,TopologyInfo> fMapInfo;
    #endif
};

#endif
