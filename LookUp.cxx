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

bool LookUp::CheckIntegrity(const AliITSMFTClusterPix& clust, const unsigned long& hash){
  fTop.SetPattern(clust);
  if(fTop.GetHash()!=hash){
    cout << "Different hashes" << endl;
  }
  auto ret = fDict.fFinalMap.find(fTop.GetHash());
  if(ret!=fDict.fFinalMap.end()) return ret->second;
  else{
    int rs = fTop.GetRowSpan();
    int cs = fTop.GetColumnSpan();
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
