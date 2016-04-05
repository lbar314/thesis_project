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
  if(gSystem->CompileMacro("MinimTopology.cxx",opt.data()))
    cout << endl << " MinimTopology: COMPILED" <<  endl << endl;
  else{
    cout << endl << " MinimTopology: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("Dictionary.cxx",opt.data()))
    cout << endl << " Dictionary: COMPILED" <<  endl << endl;
  else{
    cout << endl << " Dictionary: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("BuildDictionary.cxx",opt.data()))
    cout << endl << " BuildDictionary: COMPILED" <<  endl << endl;
  else{
    cout << endl << " BuildDictionary: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("LookUp.cxx",opt.data()))
    cout << endl << " LookUp: COMPILED" <<  endl << endl;
  else {
    cout << endl << " LookUp: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("FastSimulation.cxx",opt.data()))
    cout << endl << " FastSimulation: COMPILED" <<  endl << endl;
  else{
    cout << endl << " FastSimulation: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("testBuild.C",opt.data()))
    cout << endl << " testBuild: COMPILED" <<  endl << endl;
  else{
    cout << endl << " testBuild: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("testSimulation.C",opt.data()))
    cout << endl << " testSimulation: COMPILED" <<  endl << endl;
  else{
    cout << endl << " testSimulation: FAILED" <<  endl << endl;
    return;
  }
}
