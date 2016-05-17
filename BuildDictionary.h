#ifndef BUILDDICTIONARY_H
#define BUILDDICTIONARY_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "./MinimTopology.h"
#include "./Dictionary.h"
//#include "TH1F.h"

#define _STUDY_
//#define _HISTO_ //in order to have a histogram with the ditribution of groupIDs

#ifdef _STUDY_
  struct TopologyInfo{
    int sizeX;
    int sizeZ;
    float xCOG;
    float zCOG;
    float xMean;
    float xSigma2;
    float zMean;
    float zSigma2;
    int nPixels;
  };
#endif

struct groupTmp{
  unsigned long GrCounts;
  int tempGroupID;
  vector<unsigned long> GrHashes;
  int groupID;
};

class BuildDictionary {

  public:
    #ifdef _HISTO_
      TH1F fHdist; //Distribution of groupIDs
    #endif

    #ifndef _STUDY_
      void AccountTopology(const std::string &cluster);
    #else
      void AccountTopology(const std::string &cluster, float dX, float dZ);
    #endif

    unsigned long checkHash(const std::string& clust);

    BuildDictionary();

    void SetNGroups(unsigned int ngr); //Set number of groups
    void SetThreshold(double thr);
    void SetThresholdCumulative(double cumulative); //Considering the integral
    void Grouping();
    std::ostream& showMap(std::ostream &out);
    void PrintDictionary(string fname);

    int GetTotClusters() const {return fTotClusters;}
    int GetNotInGroups() const {return fNotInGroups;}
    int GetNGroups() const {return fNGroups;}

    Dictionary fDict;

  private:
    map<unsigned long,pair<MinimTopology,unsigned long>> fMapTop; //<hash,<topology,counts>>,
    vector <pair<unsigned long,unsigned long>> fTopFreq; //<freq,hash>, needed to define threshold
    int fTotClusters;
    int fNGroups;
    int fNotInGroups;
    double fThreshold;
    #ifdef _STUDY_
      map<long unsigned,TopologyInfo> fMapInfo;
    #endif
    MinimTopology fTop;
};

#endif
