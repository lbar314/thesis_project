#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <iostream>
#endif

void compiler(TString myopt="fast"){
  string opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("Topology.cxx",opt.data());
  cout << endl << " Topology: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("TopDatabase.cxx",opt.data());
  cout << endl << " TopDatabase: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("MinimTopology.cxx",opt.data());
  cout << endl << " MinimTopology: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("MinimDatabase.cxx",opt.data());
  cout << endl << " MinimDatabase: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("Dictionary.cxx",opt.data());
  cout << endl << " Dictionary: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("BuildDictionary.cxx",opt.data());
  cout << endl << " BuildDictionary: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("LookUp.cxx",opt.data());
  cout << endl << " LookUp: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("FastSimulation.cxx",opt.data());
  cout << endl << " FastSimulation: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("debug.C",opt.data());
  cout << endl << " debug: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("debugMin.C",opt.data());
  cout << endl << " debugMin: COMPILED" <<  endl << endl;
  gSystem->CompileMacro("testSimulation.C",opt.data());
  cout << endl << " testSimulation: COMPILED" <<  endl << endl;

}
