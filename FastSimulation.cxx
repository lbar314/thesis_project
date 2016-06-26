#include "./FastSimulation.h"
#include <random>
#include <algorithm>
#include <iostream>

FastSimulation::FastSimulation(std::string fname,unsigned seed){
  fDict.ReadBinary(fname);
  generator = std::mt19937(seed);
  distribution = std::uniform_real_distribution<double>(0.0,1.0);

}

int FastSimulation::GetRandom(){
  double rnd = distribution(generator);
  auto ind = std::upper_bound(fDict.fGroupVec.begin(),fDict.fGroupVec.end(),rnd, [] (const double &comp1, const GroupStr &comp2) {return comp1<comp2.freq;});
  return std::distance(fDict.fGroupVec.begin(),ind);
}
