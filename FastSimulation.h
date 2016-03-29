#ifndef FASTSIMULATION_H
#define FASTSIMULATION_H
#include "./Dictionary.h"

class FastSimulation{

  public:
    FastSimulation(std::string fname);
    int GetRandom();

  private:
    Dictionary fDict;
};

#endif
