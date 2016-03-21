#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "./MinimTopology.h"
#include "./MinimDatabase.h"
#include "TRandom.h"

using namespace std;

MinimDatabase::MinimDatabase():fMapTop(),fFinalMap(),fTotClusters(0),fNGroups(0),fTopFreq(),fNotInGroups(0)
#ifdef _HISTO_
  ,fHdist()
#endif
#ifdef _STUDY_
  ,fMapInfo()
#endif
{
}

MinimDatabase::~MinimDatabase(){
}

#ifndef _STUDY_
  void MinimDatabase::AccountTopology(const AliITSMFTClusterPix &cluster){
    fTotClusters++;
    MinimTopology top(cluster);
    top.GetHash();
    //pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
    auto ret = fMapTop.insert(make_pair(top.GetHash(),make_pair(top,1)));
    if(ret.second==false) ret.first->second.second++;
  }
#endif

#ifdef _STUDY_
  void MinimDatabase::AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ){
    fTotClusters++;
    MinimTopology top(cluster);
    //pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
    auto ret = fMapTop.insert(make_pair(top.GetHash(),make_pair(top,1)));
    if(ret.second==true){
      //___________________DEFINING_TOPOLOGY_CHARACTERISTICS__________________
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
      for(unsigned int i=2; i<top.GetPattern().length(); i++){
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
      topInf.xCOG = 0.5 + (float)tempyCOG/(float)tempFiredPixels;
      topInf.zCOG = 0.5 + (float)tempzCOG/(float)tempFiredPixels;
      topInf.nPixels = tempFiredPixels;
      topInf.hdX = TH1F(Form("hdX%lu",top.GetHash()),"#DeltaX",10,-5e-3,5e-3);
      topInf.hdX.Fill(dX);
      topInf.hdZ = TH1F(Form("hdA%lu",top.GetHash()),"#DeltaZ",10,-5e-3,5e-3);
      topInf.hdZ.Fill(dZ);
      fMapInfo.insert(make_pair(top.GetHash(),topInf));
    }
    else{
      ret.first->second.second++;
      auto ind = fMapInfo.find(top.GetHash());
      ind->second.hdX.Fill(dX);
      ind->second.hdZ.Fill(dZ);
    }
  }
#endif

void MinimDatabase::SetThresholdCumulative(float cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  float totFreq = 0.;
  for(auto &&p : fMapTop){
    fTopFreq.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(fTopFreq.begin(),fTopFreq.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  fNotInGroups = 0;
  fNGroups = 0;
  fFinalMap.clear();
  for(auto &q : fTopFreq){
    totFreq += (q.first)/fTotClusters;
    if(totFreq<cumulative){
      fNotInGroups++;
      fFinalMap.insert(make_pair(q.second,fNGroups++));
    }
    else break;
  }
}

void MinimDatabase::Grouping(){

  #ifdef _HISTO_
    fHdist = TH1F("fHdist", "Groups distribution", fNGroups+4, -0.5, fNGroups+3.5);
    fHdist.GetXaxis()->SetTitle("GroupID");
    fHdist.SetFillColor(kRed);
    fHdist.SetFillStyle(3005);

    for(int j=0; j<fNotInGroups; j++){
      fHdist.Fill(j,fTopFreq[j].first);
      //rough estimation fo the error considering a uniform distribution
      //fGroupVec[j].errX = (fMapTop.find(hash)->second.first.GetRowSpan())/(cmath::sqrt(12));
      //fGroupVec[j].errZ = (fMapTop.find(hash)->second.first.GetColumnSpan())/(cmath::sqrt(12));
    }
  #endif
  //This is just a dummy grouping
  for(unsigned int j = (unsigned int)fNotInGroups; j<fTopFreq.size(); j++){
    unsigned long int &hash = fTopFreq[j].second;
    int rs = fMapTop.find(hash)->second.first.GetRowSpan();
    int cs = fMapTop.find(hash)->second.first.GetColumnSpan();
    int box = rs*cs;
    if(box < 10){
      fFinalMap.insert(make_pair(hash,fNGroups));
      #ifdef _HISTO_
        fHdist.Fill(fNGroups,fTopFreq[j].first);
      #endif
    }
    else if(box >=10 && box<25){
      fFinalMap.insert(make_pair(hash,fNGroups+1));
      #ifdef _HISTO_
        fHdist.Fill(fNGroups+1,fTopFreq[j].first);
      #endif
    }
    else if(box >=25 && box<50){
      fFinalMap.insert(make_pair(hash,fNGroups+2));
      #ifdef _HISTO_
        fHdist.Fill(fNGroups+2,fTopFreq[j].first);
      #endif
    }
    else{
      fFinalMap.insert(make_pair(hash,fNGroups+3));
      #ifdef _HISTO_
        fHdist.Fill(fNGroups+3,fTopFreq[j].first);
      #endif
    }
  }
  fNGroups+=4;
}

std::ostream& MinimDatabase::showMap(std::ostream &out){
  for(auto &p : fMapTop){
    out << "Hash: " << p.second.first.GetHash() << endl;
    p.second.first.printTop(out);
  }
}
