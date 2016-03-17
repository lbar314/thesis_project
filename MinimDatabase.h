#ifndef MINIMDATABASE_H
#define MINIMDATABASE_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TH1F.h"
#include <map>
#include "AliITSMFTClusterPix.h"
#include "./MinimTopology.h"

#define _STUDY_

class MinimDatabase {

  public:
    TH1F fHdist;

    MinimDatabase();
    ~MinimDatabase();

    #ifndef _STUDY_
      void AccountTopology(const AliITSMFTClusterPix &cluster);
    #endif

    #ifdef _STUDY_
      void AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ);
    #endif

    void SetThresholdCumulative(float cumulative);
    std::ostream& showMap(std::ostream &out);

    void Grouping();
    int GetTotClusters() const {return fTotClusters;}
    int GetNotInGroups() const {return fNotInGroups;}
    int GetNGroups() const {return fNGroups;}

  private:
    map<unsigned long,pair<MinimTopology,unsigned long>> fMapTop; //<hash,<topology,counts>>,
    map<unsigned long, int> fFinalMap; //<hash,groupID>
    vector <pair<unsigned long,unsigned long>> freqv; //<freq,hash>, needed during DB construction
    int fTotClusters;
    int fNGroups;
    int fNotInGroups;

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
    map<long unsigned,TopologyInfo> fMapInfo;

  #endif
};

#endif
