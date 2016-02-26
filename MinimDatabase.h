#ifndef MINIMDATABASE_H
#define MINIMDATABASE_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "AliITSMFTClusterPix.h"
#include "./MinimTopology.h"

class MinimDatabase{

  public:
    MinimDatabase();
    virtual ~MinimDatabase();
    void AccountTopology(const AliITSMFTClusterPix &cluster);
    void SetThresholdCumulative(float cumulative);
    static bool countsCompare(const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){
      return (couple1.first > couple2.first);
    }

  private:
    map<unsigned long,pair<MinimTopology,unsigned long>> fMapTop; //<hash,<topology,counts>>,
    map<unsigned long, int> fFinalMap; //<hash,groupID>
    vector <pair<unsigned long,unsigned long>> freqv; //<freq,hash>, needed during DB construction
    int fTotClusters;
    int fNGroups;
    int fNotInGroups;
};

#endif
