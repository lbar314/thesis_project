#include "./LookUp.h"

LookUp::LookUp(std::string fname){
  fDict.ReadFile(fname);
  fOver = fDict.fFinalMap.size();
}

int LookUp::GroupFinder(const AliITSMFTClusterPix& clust){
  fTop.SetPattern(clust);
  auto ret = fDict.fFinalMap.find(fTop.GetHash());
  if(ret!=fDict.fFinalMap.end()) return ret->second;
  else{
    int rs = fTop.GetRowSpan();
    int cs = fTop.GetColumnSpan();
    int index = (rs/5)*7 + cs/5;
    if(index >48) index = 48;
    return (fOver+index);
  }
}
