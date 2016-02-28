#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include "./MinimTopology.h"
#include "./MinimDatabase.h"

using namespace std;

MinimDatabase::MinimDatabase():fMapTop(),fFinalMap(),fTotClusters(0),fNGroups(0),freqv(),fNotInGroups(0){
}

MinimDatabase::~MinimDatabase(){
}

void MinimDatabase::AccountTopology(const AliITSMFTClusterPix &cluster){
  fTotClusters++;
  MinimTopology top(cluster);
  pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
  ret = fMapTop.insert(make_pair(top.GetHash(),make_pair(top,1)));
  if(ret.second==false) ret.first->second.second++;
}

void MinimDatabase::SetThresholdCumulative(float cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  float totFreq = 0.;
  for(auto &&p : fMapTop){
    freqv.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(freqv.begin(),freqv.end(), MinimDatabase::countsCompare);
  fNotInGroups = 0;
  fNGroups = 0;
  fFinalMap.clear();
  for(auto &q : freqv){
    totFreq += (q.first)/fTotClusters;
    if(totFreq<cumulative){
      fNotInGroups++;
      fFinalMap.insert(make_pair(q.second,fNGroups++));
    }
  }
}
