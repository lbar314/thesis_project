#ifndef DICTIONARY_H
#define DICTIONARY_H
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

using std::vector;
using std::map;
using std::string;

struct GroupStr{
  unsigned long hash;
  float errX;
  float errZ;
  double freq;
};

struct Dictionary{
  map<unsigned long, int> fFinalMap; //<hash,groupID> just for topologies over threshold
  vector<GroupStr> fGroupVec;
  friend std::ostream& operator<<(std::ostream& os, const Dictionary& dict);
  void ReadFile(string fname);
};

#endif
