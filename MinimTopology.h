#ifndef MINIMTOPOLOGY_H
#define MINIMTOPOLOGY_H
#include "AliITSMFTClusterPix.h"
#include <iostream>
#include <string>

using namespace std;

class MinimTopology {

  public:
    MinimTopology();
    MinimTopology(const AliITSMFTClusterPix &cluster);
    
    string& GetPattern() {return fPattern;}
    int GetRowSpan() const {return fPattern[0];}
    int GetColumnSpan() const {return fPattern[1];}
    unsigned long GetHash() const {return fHash;}

    std::ostream& printTop(std::ostream &out);
    static unsigned int hashFunction(const void * key, int len);

  private:
    string fPattern;
    unsigned long fHash;

};

#endif
