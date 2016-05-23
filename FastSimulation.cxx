#include "./FastSimulation.h"
#include <random>
#include <algorithm>

FastSimulation::FastSimulation(std::string fname){
  fDict.ReadBinary(fname);
}

int FastSimulation::GetRandom(){
  double &max1 = fDict.fGroupVec.back().freq;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double rnd = distribution(generator);
  auto ind = std::upper_bound(fDict.fGroupVec.begin(),fDict.fGroupVec.end(),rnd, [] (const double &comp1, const GroupStr &comp2) {return comp1<comp2.freq;});
  return std::distance(fDict.fGroupVec.begin(),ind);
}
