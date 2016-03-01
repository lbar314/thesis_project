#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "./MinimTopology.h"
#include "./MinimDatabase.h"

using namespace std;

MinimDatabase::MinimDatabase():fMapTop(),fFinalMap(),fTotClusters(0),fNGroups(0),freqv(),fNotInGroups(0)
#ifdef _STUDY_
  ,fMapInfo()
#endif
{
}

MinimDatabase::~MinimDatabase(){
}

void MinimDatabase::AccountTopology(const AliITSMFTClusterPix &cluster/*, ostream& output*/){
  fTotClusters++;
  MinimTopology top(cluster);
  top.GetHash();
  //pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
  auto ret = fMapTop.insert(make_pair(top.GetHash(),make_pair(top,1)));
  if(ret.second==false) ret.first->second.second++;

  #ifdef _STUDY_
    TopologyInfo topInf;
    int &rs = topInf.sizeX = top.GetRowSpan();
    int &cs = topInf.sizeZ = top.GetColumnSpan();
    //__________________COG_Deterrmination_____________
    int tempyCOG = 0;
    int tempzCOG = 0;
    int tempFiredPixels = 0;
    unsigned char tempChar = 0;
  	int s = 0;
  	int ic = 0;
    int ir = 0;
    for(int i=2; i<top.GetPattern().length(); i++){
  		tempChar = top.GetPattern()[i];
  		s=128;//0b10000000
      while(s>0){
        if((tempChar&s)!=0){
          tempFiredPixels++;
				  tempyCOG+=ir;
          tempzCOG+=ic;
        }
        ic++;
        s/=2;
        if((ir+1)*ic==(rs*cs)) break;
  			if(ic==cs){
          ic=0;
          ir++;
        }
      }
  		if((ir+1)*ic==(rs*cs)) break;
    }
    topInf.Y = 0.5 + (float)tempyCOG/(float)tempFiredPixels;
    topInf.Z = 0.5 + (float)tempzCOG/(float)tempFiredPixels;
    topInf.nPixels = tempFiredPixels;

    //__________________SETTING_ERRORS_________________
    topInf.sigmaY2 = cluster.GetSigmaY2();
    topInf.sigmaZ2 = cluster.GetSigmaZ2();
    topInf.sigmaYZ = cluster.GetSigmaYZ();
    topInf.sizeY = top.GetRowSpan();
    topInf.sizeZ = top.GetColumnSpan();
    pair<map<unsigned long, TopologyInfo>::iterator,bool> ind;
    ind = fMapInfo.insert(make_pair(top.GetHash(),topInf));
    if(ind.second==true){top.printTop(output);}
    // else{
    //   if(abs(topInf.sigmaY2-ind.first->second.sigmaY2)>abs(topInf.sigmaY2)*1e-7){
    //     cout << "Different sigmaY2: " << topInf.sigmaY2 - ind.first->second.sigmaY2 << endl;
    //     top.printTop(cout);
    //     cin.get();
    //   }
    //   if(abs(topInf.sigmaZ2-ind.first->second.sigmaZ2)>abs(topInf.sigmaZ2)*1e-7){
    //     cout << "Different sigmaZ2: " << topInf.sigmaZ2 - ind.first->second.sigmaZ2 << endl;
    //     top.printTop(cout);
    //     cin.get();
    //   }
    //   if(abs(topInf.sigmaYZ-ind.first->second.sigmaYZ)>abs(topInf.sigmaYZ)*1e-7){
    //     cout << "Different sigmaYZ: " << topInf.sigmaYZ - ind.first->second.sigmaYZ << endl;
    //     top.printTop(cout);
    //     cin.get();
    //   }
    // }
  #endif
}

void MinimDatabase::SetThresholdCumulative(float cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  float totFreq = 0.;
  for(auto &&p : fMapTop){
    freqv.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(freqv.begin(),freqv.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  fNotInGroups = 0;
  fNGroups = 0;
  fFinalMap.clear();
  for(auto &q : freqv){
    totFreq += (q.first)/fTotClusters;
    if(totFreq<cumulative){
      fNotInGroups++;
      fFinalMap.insert(make_pair(q.second,fNGroups++));
    }
  }
}

std::ostream& MinimDatabase::showMap(std::ostream &out){
  for(auto &p : fMapTop){
    out << "Hash: " << p.second.first.GetHash() << endl;
    p.second.first.printTop(out);
  }
}
