#include "Dictionary.h"

std::ostream& operator<<(std::ostream& os, const Dictionary& dict)
{
  for(auto &p : dict.fGroupVec){
    os << p.hash << " " << p.errX << " " << p.errZ << " " << p.freq << std::endl;
  }
  return os;
}

void Dictionary::ReadFile(string fname){
  fGroupVec.clear();
  fFinalMap.clear();
  std::ifstream in(fname);
  GroupStr gr;
  int groupID=0;
  if(!in.is_open()){
    std::cout << "The file could not be opened" << endl;
    exit(1);
  }
  else{
    while(in >> gr.hash >> gr.errX >> gr.errZ >> gr.freq){
      fGroupVec.push_back(gr);
      if(((gr.hash)&0xffffffff) != 0) fFinalMap.insert(make_pair(gr.hash,groupID));
      groupID++;
    }
  }
  in.close();
}
