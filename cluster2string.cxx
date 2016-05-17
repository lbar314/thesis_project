#include "./cluster2string.h"
#include "AliITSMFTClusterPix.h"
#include <string>

void FromCluster2String(const AliITSMFTClusterPix &cluster, std::string &str){
  int rs = cluster.GetPatternRowSpan();
  int cs = cluster.GetPatternColSpan();
  int nBytes = (rs*cs)>>3;
  if(((rs*cs)%8)!=0) nBytes++;
  str.resize(nBytes+2,0);
  str[0]=rs;
	str[1]=cs;
  cluster.GetPattern(&str[2],nBytes);
}
