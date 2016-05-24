#ifndef LOOKUP_H
#define LOOKUP_H
#include "./Dictionary.h"
#include "./MinimTopology.h"
#include <iostream>
#include <unordered_map>
#include <vector>

using std::vector;
using std::unordered_map;

class LookUp{
  public:
    LookUp(std::string fname);
    int GroupFinder(const std::string& clust);
    int GetOver() {return fOver;}

  private:
    Dictionary fDict;
    int fOver;
    MinimTopology fTop;
};

#endif
