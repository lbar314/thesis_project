#include "./cluster2string.h"
//#include "AliITSMFTClusterPix.h"
//#include "ITSMFTReconstruction/Cluster.h"

//#include <string>

void FromCluster2String(const Cluster &cluster, std::string &str){
  int rs = cluster.getPatternRowSpan();
  int cs = cluster.getPatternColSpan();
  int nBytes = (rs*cs)>>3;
  if(((rs*cs)%8)!=0) nBytes++;
  str.resize(nBytes+2,0);
  str[0]=rs;
	str[1]=cs;
  cluster.getPattern(&str[2],nBytes);
}
