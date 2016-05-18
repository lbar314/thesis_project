#include "Dictionary.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

std::ostream& operator<<(std::ostream& os, const Dictionary& dict)
{
  for(auto &p : dict.fGroupVec){
    os << p.hash << " " << p.errX << " " << p.errZ << " " << p.freq << std::endl;
  }
  return os;
}

void Dictionary::WriteBinary(string outputfile){
  std::ofstream file_output(outputfile, std::ios::out | std::ios::binary);
  for(auto &p : fGroupVec){
    file_output.write(reinterpret_cast<char *>(&p.hash),sizeof(unsigned long));
    file_output.write(reinterpret_cast<char *>(&p.errX),sizeof(float));
    file_output.write(reinterpret_cast<char *>(&p.errZ),sizeof(float));
    file_output.write(reinterpret_cast<char *>(&p.freq),sizeof(double));
  }
  file_output.close();
}

void Dictionary::ReadFile(string fname){
  fGroupVec.clear();
  fFinalMap.clear();
  std::ifstream in(fname);
  GroupStr gr;
  int groupID=0;
  if(!in.is_open()){
    cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in >> gr.hash >> gr.errX >> gr.errZ >> gr.freq){
      fGroupVec.push_back(gr);
      if(((gr.hash)&0xffffffff) != 0) fFinalMap.insert(std::make_pair(gr.hash,groupID));
      groupID++;
    }
  }
  in.close();
}

void Dictionary::ReadBinary(string fname){
  fGroupVec.clear();
  fFinalMap.clear();
  std::ifstream in(fname.data(),std::ios::in | std::ios::binary);
  GroupStr gr;
  int groupID=0;
  if(!in.is_open()){
    cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in.read(reinterpret_cast<char*>(&gr.hash), sizeof(unsigned long))){
      in.read(reinterpret_cast<char*>(&gr.errX), sizeof(float));
      in.read(reinterpret_cast<char*>(&gr.errZ), sizeof(float));
      in.read(reinterpret_cast<char*>(&gr.freq), sizeof(double));
      fGroupVec.push_back(gr);
      if(((gr.hash)&0xffffffff) != 0) fFinalMap.insert(std::make_pair(gr.hash,groupID));
      groupID++;
    }
  }
  in.close();
}
