#include "./FastSimulation.h"
#include "TRandom3.h"
#include <algorithm>

FastSimulation::FastSimulation(std::string fname){
  fDict.ReadFile(fname);
}

int FastSimulation::GetRandom(){
  const double ciccio = fDict.fGroupVec.back().freq;
  double rnd = gRandom->Uniform(ciccio);
  auto ind = std::upper_bound(fDict.fGroupVec.begin(),fDict.fGroupVec.end(),rnd, [] (const double &comp1, const GroupStr &comp2) {return comp1<comp2.freq;});
  return std::distance(fDict.fGroupVec.begin(),ind);
}
