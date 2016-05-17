#include "./MinimTopology.h"
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

vector<string> vec;

void read_clusters(string inputfile);

void timer(string dictionary_file);

void group_storer(string dictionary_file, string outut_file);

int main(){
  read_clusters("clusterlist.bin");
  timer("dizionario.txt");
  group_storer("dizionario.txt", "groups.txt");
  return 0;
}

void read_clusters(string inputfile){
  cout << "Reading clusters from file ..." << endl;
  string str;
  std::ifstream in(inputfile.data(),std::ios::in | std::ios::binary);
  std::ofstream check("bench.txt");
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
      check << ntotclusters << ") string size: " << size1 << endl;
      buf = new char[size1];
	    in.read(buf, size1);
      str.append(buf, size1);
      MinimTopology a(str);
      a.printTop(check);
      vec.push_back(str);
      delete[] buf;
      ntotclusters++;
    }
  }
  in.close();
  check.close();
  cout << "Done! Number of clusters: " << vec.size() << endl;
}

void timer(string dictionary_file){
  LookUp finder(dictionary_file);
  string str;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  for(auto &p : vec){
    str = p;
    finder.GroupFinder(p);
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
  cout << "Time for processing all the clusters: " << time_span.count() << " s" << endl;
}

void group_storer(string dictionary_file, string outut_file){
  cout << "Writing groupID to file ...";
  LookUp finder(dictionary_file);
  std::ofstream out(outut_file);
  string str;
  for(auto &p : vec){
    str = p;
    out << finder.GroupFinder(p) << endl;
  }
  out.close();
  cout << " Done!" << endl;
}
