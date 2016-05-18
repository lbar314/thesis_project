#include "./FastSimulation.h"
#include "TRandom3.h"
#include <algorithm>

FastSimulation::FastSimulation(std::string fname){
  fDict.ReadBinary(fname);
}

int FastSimulation::GetRandom(){
  double &max1 = fDict.fGroupVec.back().freq;
  double rnd = gRandom->Uniform(max1);
  auto ind = std::upper_bound(fDict.fGroupVec.begin(),fDict.fGroupVec.end(),rnd, [] (const double &comp1, const GroupStr &comp2) {return comp1<comp2.freq;});
  return std::distance(fDict.fGroupVec.begin(),ind);
}
