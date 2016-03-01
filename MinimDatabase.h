#ifndef MINIMDATABASE_H
#define MINIMDATABASE_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "AliITSMFTClusterPix.h"
#include "./MinimTopology.h"

//#define _STUDY_

class MinimDatabase {

  public:
    MinimDatabase();
    ~MinimDatabase();
    void AccountTopology(const AliITSMFTClusterPix &cluster/*, ostream& output=cout*/);
    void SetThresholdCumulative(float cumulative);
    std::ostream& showMap(std::ostream &out);

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
      float X;
      float Z;
      float sigmaY2;
      float sigmaZ2;
      float sigmaYZ;
      int nPixels;
    };
    map<long unsigned,TopologyInfo> fMapInfo;

  #endif
};

#endif
