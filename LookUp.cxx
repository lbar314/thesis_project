#include "./LookUp.h"

LookUp::LookUp(std::string fname){
  fDict.ReadFile(fname);
  fOver = fDict.fFinalMap.size();
}

int LookUp::GroupFinder(AliITSMFTClusterPix& clust){
  MinimTopology top(clust);
  auto ret = fDict.fFinalMap.find(top.GetHash());
  if(ret!=fDict.fFinalMap.end()) return ret->second;
  else{
    int rs = top.GetRowSpan();
    int cs = top.GetColumnSpan();
    int box = rs*cs;
    if(box < 10){
      return fOver;
    }
    else if(box >=10 && box<25){
      return fOver+1;
    }
    else if(box >=25 && box<100){
      return fOver+2;
    }
    else{
     return fOver+3;
    }
  }
}
