#ifndef MINIMTOPOLOGY_H
#define MINIMTOPOLOGY_H
#include <iostream>
#include <string>

class MinimTopology {

  public:
    MinimTopology();
    MinimTopology(const std::string &str);

    std::string& GetPattern() {return fPattern;}
    int GetRowSpan() const {return fPattern[0];}
    int GetColumnSpan() const {return fPattern[1];}
    unsigned long GetHash() const {return fHash;}
    friend std::ostream& operator<<(std::ostream& os, const MinimTopology& top);
    static unsigned int hashFunction(const void * key, int len);
    void SetPattern(const std::string &str);

  private:

    std::string fPattern;
    unsigned long fHash;

};

#endif
