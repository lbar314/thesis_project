#ifndef LOOKUP_H
#define LOOKUP_H
#include "./Dictionary.h"
#include "./MinimTopology.h"
#include "AliITSUClusterPix.h"
#include <iostream>
#include <map>
#include <vector>

using std::vector;
using std::map;

class LookUp{
  public:
    LookUp(std::string fname);
    int GroupFinder(AliITSMFTClusterPix& clust);

  private:
    Dictionary fDict;
    int fOver;

};

#endif
