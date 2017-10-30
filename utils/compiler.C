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
  if(gSystem->CompileMacro("cluster2string.cxx",opt.data()))
    cout << endl << " cluster2string: COMPILED" <<  endl << endl;
  else{
    cout << endl << " cluster2string: FAILED" <<  endl << endl;
    return;
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
  if(gSystem->CompileMacro("testLookUp.C",opt.data()))
    cout << endl << " testLookUp: COMPILED" <<  endl << endl;
  else{
    cout << endl << " testLookUp: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("StoreClusters.C",opt.data()))
    cout << endl << " StoreClusters: COMPILED" <<  endl << endl;
  else{
    cout << endl << " StoreClusters: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("timeDistr.C",opt.data()))
    cout << endl << " timeDistr: COMPILED" <<  endl << endl;
  else{
    cout << endl << " timeDistr: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("CheckChain.C",opt.data()))
    cout << endl << " CheckChain: COMPILED" <<  endl << endl;
  else{
    cout << endl << " CheckChain: FAILED" <<  endl << endl;
    return;
  }
  if(gSystem->CompileMacro("HistoFiller.C",opt.data()))
    cout << endl << " HistoFiller: COMPILED" <<  endl << endl;
  else{
    cout << endl << " HistoFiller: FAILED" <<  endl << endl;
    return;
  }
}
