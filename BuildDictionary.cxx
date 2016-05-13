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

  BuildDictionary::BuildDictionary():fTotClusters(0){}

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

unsigned long BuildDictionary::checkHash(const AliITSMFTClusterPix& clust){
  fTop.SetPattern(clust);
  return fTop.GetHash();
}

void BuildDictionary::SetThreshold(double thr){
  for(auto &&p : fMapTop){
    fTopFreq.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(fTopFreq.begin(),fTopFreq.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  fNotInGroups = 0;
  fNGroups = 0;
  fDict.fFinalMap.clear();
  fThreshold=thr;
  for(auto &q : fTopFreq){
    if( ((double)q.first)/fTotClusters > thr ) fNotInGroups++;
    else break;
  }
  fNGroups=fNotInGroups;
}

void BuildDictionary::SetNGroups(unsigned int ngr){
  for(auto &&p : fMapTop){
    fTopFreq.push_back(make_pair(p.second.second,p.first));
  }
  std::sort(fTopFreq.begin(),fTopFreq.end(), [] (const pair<unsigned long, unsigned long> &couple1, const pair<unsigned long, unsigned long> &couple2){return (couple1.first > couple2.first);});
  if(ngr<10 || ngr > (fTopFreq.size()-49)){
    cout << "BuildDictionary::SetNGroups : Invalid number of groups" << endl;
    exit(1);
  }
  fNGroups = fNotInGroups = ngr-49;
  fDict.fFinalMap.clear();
  fThreshold=((double)fTopFreq[fNotInGroups-1].first)/fTotClusters;
}

void BuildDictionary::SetThresholdCumulative(double cumulative){
  cout<<"SetThresholdCumulative: fTotClusters: " << fTotClusters << endl;
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

  cout<<"Grouping: fTotClusters: " << fTotClusters << endl;
  #ifdef _HISTO_
    fHdist = TH1F("fHdist", "Groups distribution", fNGroups+49, -0.5, fNGroups+48.5);
    fHdist.GetXaxis()->SetTitle("GroupID");
    fHdist.SetFillColor(kRed);
    fHdist.SetFillStyle(3005);
  #endif

  double totFreq=0.;
  for(int j=0; j<fNotInGroups; j++){
    #ifdef _HISTO_
      fHdist.Fill(j,fTopFreq[j].first);
    #endif
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
  //Grouping based on binning over number of rows and columns (7*7)
  fNGroups+=49; //(7*7)
  //array of groups
  std::array<GroupStr,49> GroupArray;
  std::array<unsigned long,49> groupCounts{0};
  auto func = [&GroupArray] (int rowBinEdge, int colBinEdge, int &index) {
    unsigned long provvHash = 0;
    provvHash = ( ((unsigned long)(index+1)) << 32 ) & 0xffffffff00000000;
    GroupArray[index].hash = provvHash;
    GroupArray[index].errX = (rowBinEdge)*2e-3/std::sqrt(12); // 2e-3 is the pitch
    GroupArray[index].errZ = (colBinEdge)*2e-3/std::sqrt(12); // 2e-3 is the pitch
    index++;
    return;
  };
  int grNum=0;
  for(int ir=0; ir<6; ir++){ //row bins: {[0;4],[5;9],[10;14],[15;19],[20;24],[25,29]} (+ [30;32] later)
    for(int ic=0; ic<6; ic++){ //col bins: {[0;4],[5;9],[10;14],[15;19],[20;24],[25,29]} (+ [30;32] later)
      func((ir+1)*5-1, (ic+1)*5-1, grNum);
    }
    // col bin [30;32]
    func((ir+1)*5-1, 32, grNum);
  }
  // row bin [30;32]
  for(int ic=0; ic<6; ic++){ //col bins: {[0;4],[5;9],[10;14],[15;19],[20;24],[25,29]} (+ [30;32] later)
    func(32, (ic+1)*5-1, grNum);
    unsigned long provvHash = 0;
  }
  func(32, 32, grNum);
  if(grNum!=49){
    cout << "Wrong number of groups" << endl;
    exit(1);
  }

  cout << endl;unsigned long hash1;
  int rs;
  int cs;
  int index;

  for(unsigned int j = (unsigned int)fNotInGroups; j<fTopFreq.size(); j++){
    unsigned long
    hash1 = fTopFreq[j].second;
    rs = fMapTop.find(hash1)->second.first.GetRowSpan();
    cs = fMapTop.find(hash1)->second.first.GetColumnSpan();
    index = (rs/5)*7 + cs/5;
    if(index >48) index = 48;
    groupCounts[index]+=fTopFreq[j].first;
  }

  for(int i=0; i<49; i++){
    totFreq+=((double)groupCounts[i])/fTotClusters;
    GroupArray[i].freq = totFreq;
    #ifdef _HISTO_
      fHdist.Fill(fNotInGroups+i,groupCounts[i]);
    #endif
    fDict.fGroupVec.push_back(GroupArray[i]);
  }
  #ifdef _HISTO_
    fHdist.Scale(1./fHdist.Integral());
  #endif
}

std::ostream& BuildDictionary::showMap(std::ostream &out){
  out << "Vecchia versione" << endl;
  for(auto &p : fMapTop){
    out << "Hash: " << p.second.first.GetHash() << endl;
    out << "counts: " << p.second.second << endl;
    p.second.first.printTop(out);
  }
}

void BuildDictionary::PrintDictionary(string fname){
  ofstream out(fname);
  out << fDict;
  out.close();
}
