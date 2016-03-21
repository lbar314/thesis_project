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
  if(cumulative<=0. || cumulative >=1.) cumulative = 0.99;
  float totFreq = 0.;
  for(auto &&p : fMapTop){
    fTopFreq.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(fTopFreq.begin(),fTopFreq.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  fNotInGroups = 0;
  fNGroups = 0;
  fFinalMap.clear();
  for(auto &q : fTopFreq){
    totFreq += ((float)(q.first))/fTotClusters;
    if(totFreq<cumulative){
      fNotInGroups++;
      //fFinalMap.insert(make_pair(q.second,fNGroups++));
    }
    else break;
  }
  fThreshold=((float)(fTopFreq[--fNotInGroups].first))/fTotClusters;

  while(((float)fTopFreq[fNotInGroups].first)/fTotClusters == fThreshold) fNotInGroups--;
  fThreshold=((float)fTopFreq[fNotInGroups++].first)/fTotClusters;
  fNGroups=fNotInGroups;
}

void MinimDatabase::Grouping(){
  struct groupTmp{
    unsigned long GrCounts;
    int tempGroupID;
    vector<unsigned long> GrHashes;
    //float tempErrX;
    //float tempErrZ;
    int groupID;
  };
  vector<groupTmp> tempGroupVector;
  for(int j=0; j<fNotInGroups; j++){
    groupTmp gr;
    gr.GrCounts = fTopFreq[j].first;
    gr.GrHashes.push_back(fTopFreq[j].second);
    gr.tempGroupID = j;
    //rough estimation fo the error considering a uniform distribution
    //fGroupVec[j].errX = (fMapTop.find(hash)->second.first.GetRowSpan())/(cmath::sqrt(12));
    //fGroupVec[j].errZ = (fMapTop.find(hash)->second.first.GetColumnSpan())/(cmath::sqrt(12));
    tempGroupVector.push_back(gr);
  }
  //This is just a dummy grouping
  fNGroups+=4; //for choice
  //group 1:
  groupTmp gr1;
    gr1.GrCounts = 0;
    gr1.tempGroupID = fNGroups;
  //group 2:
  groupTmp gr2;
    gr2.GrCounts = 0;
    gr2.tempGroupID = fNGroups+1;
  //group 3:
  groupTmp gr3;
    gr3.GrCounts = 0;
    gr3.tempGroupID = fNGroups+2;
  //group 3:
  groupTmp gr4;
    gr4.GrCounts = 0;
    gr4.tempGroupID = fNGroups+3;

  for(unsigned int j = (unsigned int)fNotInGroups; j<fTopFreq.size(); j++){
    unsigned long int &hash = fTopFreq[j].second;
    int rs = fMapTop.find(hash)->second.first.GetRowSpan();
    int cs = fMapTop.find(hash)->second.first.GetColumnSpan();
    int box = rs*cs;
    if(box < 10){
      gr1.GrCounts+=fTopFreq[j].first;
      gr1.GrHashes.push_back(hash);
    }
    else if(box >=10 && box<25){
      gr2.GrCounts+=fTopFreq[j].first;
      gr2.GrHashes.push_back(hash);
    }
    else if(box >=25 && box<50){
      gr1.GrCounts+=fTopFreq[j].first;
      gr1.GrHashes.push_back(hash);
    }
    else{
      gr1.GrCounts+=fTopFreq[j].first;
      gr1.GrHashes.push_back(hash);
    }
  }
  tempGroupVector.push_back(gr1);
  tempGroupVector.push_back(gr2);
  tempGroupVector.push_back(gr3);
  tempGroupVector.push_back(gr4);
  //sorting the temporary array
  std::sort(tempGroupVector.begin(),tempGroupVector.end(), [] (const groupTmp &comp1, const groupTmp &comp2){return (comp1.GrCounts > comp2.GrCounts);});

  #ifdef _HISTO_
    fHdist = TH1F("fHdist", "Groups distribution", fNGroups+4, -0.5, fNGroups+3.5);
    fHdist.GetXaxis()->SetTitle("GroupID");
    fHdist.SetFillColor(kRed);
    fHdist.SetFillStyle(3005);
  #endif

  int iDef=0;
  for(auto &p : tempGroupVector){
    p.groupID = iDef++;
    fHdist.Fill(p.groupID,p.GrCounts);
    for(auto &q : p.GrHashes){
      fFinalMap.insert(make_pair(q,p.groupID));
    }
  }
}

std::ostream& MinimDatabase::showMap(std::ostream &out){
  for(auto &p : fMapTop){
    out << "Hash: " << p.second.first.GetHash() << endl;
    p.second.first.printTop(out);
  }
}
