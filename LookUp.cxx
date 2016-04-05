#include "./LookUp.h"

LookUp::LookUp(std::string fname){
  fDict.ReadFile(fname);
}

int LookUp::GroupFinder(const AliITSMFTClusterPix& clust){
  MinimTopology top(clust);
  int over = fDict.fFinalMap.size();
  auto ret = fDict.fFinalMap.find(top.GetHash());
  if(ret!=fDict.fFinalMap.end()) return ret->second;
  else{
    int rs = top.GetRowSpan();
    int cs = top.GetColumnSpan();
    int box = rs*cs;
    if(box < 10){
      return over;
    }
    else if(box >=10 && box<25){
      return over+1;
    }
    else if(box >=25 && box<100){
      return over+2;
    }
    else{
     return over+3;
    }
  }
}
