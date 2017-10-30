#include "./MinimTopology.h"
#include "./BuildDictionary.h"
#include "./Dictionary.h"
#include "./LookUp.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <ratio>
#include <chrono>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace std::chrono;

vector<string> vecTop;
vector<float> vecX;
vector<float> vecZ;

void read_clusters(string inputfile);
void builder();
void group_storer(string dictionary_file, string outut_file);

int main(){
  read_clusters("CompleteList.bin");
  builder();
  group_storer("dizionario.bin", "groupsBuild.txt");
  return 0;
}

void read_clusters(string inputfile){
  cout << "Reading clusters from file ..." << endl;
  string str;
  std::ifstream in(inputfile.data(),std::ios::in | std::ios::binary);
  int ntotclusters = 0;
  char* buf;
  if(!in.is_open()){
    cout << "The file could not be opened" << endl;
    return;
  }
  else{
    int size1;
    while(in.read(reinterpret_cast<char*>(&size1), sizeof(int))){
      string str="";
      buf = new char[size1];
	    in.read(buf, size1);
      str.append(buf, size1);
      vecTop.push_back(str);
      delete[] buf;
      float dX;
      in.read(reinterpret_cast<char*>(&dX), sizeof(float));
      vecX.push_back(dX);
      float dZ;
      in.read(reinterpret_cast<char*>(&dZ), sizeof(float));
      vecZ.push_back(dZ);
      ntotclusters++;
    }
  }
  in.close();
  cout << "Done! Number of clusters: " << vecTop.size() << endl;
}

void group_storer(string dictionary_file, string outut_file){
  cout << "Writing groupID to file ...";
  LookUp finder(dictionary_file);
  std::ofstream out(outut_file);
  string str;
  for(auto &p : vecTop){
    str = p;
    out << finder.GroupFinder(p) << endl;
  }
  out.close();
  cout << " Done!" << endl;
}

void builder(){
  BuildDictionary minDB;
  for(unsigned int i=0; i<vecTop.size(); i++){
    minDB.AccountTopology(vecTop[i],vecX[i], vecZ[i]);
  }
  minDB.SetThreshold(0.00001);
  minDB.Grouping();
  minDB.PrintDictionary("dizionario.txt");
  minDB.PrintDictionaryBin("dizionario.bin");
}
