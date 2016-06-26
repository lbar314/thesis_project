#ifndef FASTSIMULATION_H
#define FASTSIMULATION_H
#include "./Dictionary.h"
#include <random>

class FastSimulation{

  public:
    FastSimulation(std::string fname, unsigned seed=0xdeadbeef);
    int GetRandom();

  private:
    Dictionary fDict;
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
};

#endif
