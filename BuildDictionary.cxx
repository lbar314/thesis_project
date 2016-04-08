#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "./MinimTopology.h"
#include "./BuildDictionary.h"
#include "TRandom.h"

using namespace std;

#ifndef _STUDY_
  void BuildDictionary::AccountTopology(const AliITSMFTClusterPix &cluster){
    fTotClusters++;
    fTtop.SetPattern(cluster);
    //pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
    auto ret = fMapTop.insert(make_pair(fTop.GetHash(),make_pair(fTop,1)));
    if(ret.second==false) ret.first->second.second++;
  }
#else
  void BuildDictionary::AccountTopology(const AliITSMFTClusterPix &cluster, float dX, float dZ){
    fTotClusters++;
    fTop.SetPattern(cluster);
    //pair<map<unsigned long, pair<MinimTopology,unsigned long>>::iterator,bool> ret;
    auto ret = fMapTop.insert(make_pair(fTop.GetHash(),make_pair(fTop,1)));
    if(ret.second==true){
      //___________________DEFINING_TOPOLOGY_CHARACTERISTICS__________________
      TopologyInfo topInf;
      int &rs = topInf.sizeX = fTop.GetRowSpan();
      int &cs = topInf.sizeZ = fTop.GetColumnSpan();
      //__________________COG_Deterrmination_____________
      int tempyCOG = 0;
      int tempzCOG = 0;
      int tempFiredPixels = 0;
      unsigned char tempChar = 0;
    	int s = 0;
    	int ic = 0;
      int ir = 0;
      for(unsigned int i=2; i<fTop.GetPattern().length(); i++){
    		tempChar = fTop.GetPattern()[i];
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
      topInf.hdX = TH1F(Form("hdX%lu",fTop.GetHash()),"#DeltaX",10,-5e-3,5e-3);
      topInf.hdX.Fill(dX);
      topInf.hdZ = TH1F(Form("hdA%lu",fTop.GetHash()),"#DeltaZ",10,-5e-3,5e-3);
      topInf.hdZ.Fill(dZ);
      fMapInfo.insert(make_pair(fTop.GetHash(),topInf));
    }
    else{
      ret.first->second.second++;
      auto ind = fMapInfo.find(fTop.GetHash());
      ind->second.hdX.Fill(dX);
      ind->second.hdZ.Fill(dZ);
    }
  }
#endif

void BuildDictionary::SetThresholdCumulative(double cumulative){
  if(cumulative<=0. || cumulative >=1.) cumulative = 0.99;
  double totFreq = 0.;
  for(auto &&p : fMapTop){
    fTopFreq.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(fTopFreq.begin(),fTopFreq.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  fNotInGroups = 0;
  fNGroups = 0;
  fDict.fFinalMap.clear();
  for(auto &q : fTopFreq){
    totFreq += ((double)(q.first))/fTotClusters;
    if(totFreq<cumulative){
      fNotInGroups++;
    }
    else break;
  }
  fThreshold=((double)(fTopFreq[--fNotInGroups].first))/fTotClusters;
  while(((double)fTopFreq[fNotInGroups].first)/fTotClusters == fThreshold) fNotInGroups--;
  fThreshold=((double)fTopFreq[fNotInGroups++].first)/fTotClusters;
  fNGroups=fNotInGroups;
}

void BuildDictionary::Grouping(){

  #ifdef _HISTO_
    fHdist = TH1F("fHdist", "Groups distribution", fNGroups+4, -0.5, fNGroups+3.5);
    fHdist.GetXaxis()->SetTitle("GroupID");
    fHdist.SetFillColor(kRed);
    fHdist.SetFillStyle(3005);
  #endif

  double totFreq=0.;
  for(int j=0; j<fNotInGroups; j++){
    fHdist.Fill(j,fTopFreq[j].first);
    totFreq+=((double)(fTopFreq[j].first))/fTotClusters;
    GroupStr gr;
    gr.hash=fTopFreq[j].second;
    gr.freq=totFreq;
    //rough estimation fo the error considering a uniform distribution
    gr.errX = (fMapTop.find(gr.hash)->second.first.GetRowSpan())/(std::sqrt(12));
    gr.errZ = (fMapTop.find(gr.hash)->second.first.GetColumnSpan())/(std::sqrt(12));
    fDict.fGroupVec.push_back(gr);
    fDict.fFinalMap.insert(make_pair(gr.hash,j));
  }
  //This is just a dummy grouping
  fNGroups+=4; //for choice
  //group 1:
  GroupStr gr1;
    gr1.hash = ((unsigned long)1) << 56;
    gr1.errX = gr1.errZ = std::sqrt(10)*2e-3/std::sqrt(12);
    unsigned long count1 = 0;
  //group 2:
  GroupStr gr2;
    gr2.hash = ((unsigned long)2) << 56;
    gr2.errX = gr2.errZ = 5.*2e-3/std::sqrt(12);
    unsigned long count2 = 0;
  //group 3:
  GroupStr gr3;
    gr3.hash = ((unsigned long)3) << 56;
    gr3.errX = gr1.errZ = 10.*2e-3/std::sqrt(12);
    unsigned long count3 = 0;
  //group 3:
  GroupStr gr4;
    gr4.hash = ((unsigned long)4) << 56;
    gr4.errX = gr1.errZ = 30.*2e-3/std::sqrt(12);
    unsigned long count4 = 0;

  for(unsigned int j = (unsigned int)fNotInGroups; j<fTopFreq.size(); j++){
    unsigned long int &hash = fTopFreq[j].second;
    int rs = fMapTop.find(hash)->second.first.GetRowSpan();
    int cs = fMapTop.find(hash)->second.first.GetColumnSpan();
    int box = rs*cs;
    if(box < 10){
      count1+=fTopFreq[j].first;
    }
    else if(box >=10 && box<25){
      count2+=fTopFreq[j].first;
    }
    else if(box >=25 && box<100){
      count3+=fTopFreq[j].first;
    }
    else{
      count4+=fTopFreq[j].first;
    }
  }
  totFreq+=((double)count1)/fTotClusters;
  gr1.freq=totFreq;
  fHdist.Fill(fNotInGroups,count1);
  fDict.fGroupVec.push_back(gr1);
  totFreq+=((double)count2)/fTotClusters;
  gr2.freq=totFreq;
  fHdist.Fill(fNotInGroups+1,count2);
  fDict.fGroupVec.push_back(gr2);
  totFreq+=((double)count3)/fTotClusters;
  gr3.freq=totFreq;
  fHdist.Fill(fNotInGroups+2,count3);
  fDict.fGroupVec.push_back(gr3);
  totFreq+=((double)count4)/fTotClusters;
  gr4.freq=totFreq;
  fHdist.Fill(fNotInGroups+3,count4);
  fDict.fGroupVec.push_back(gr4);
}

std::ostream& BuildDictionary::showMap(std::ostream &out){
  out << "Nuova versione"<< endl;
  for(auto &p : fMapTop){
    out << "Hash: " << p.second.first.GetHash() << endl;
    p.second.first.printTop(out);
  }
}

void BuildDictionary::PrintDictionary(string fname){
  ofstream out(fname);
  out << fDict;
  out.close();
}
