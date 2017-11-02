#ifndef CLUSTER2STRING_H
#define CLUSTER2STRING_H

//#include "AliITSMFTClusterPix.h"
#include "ITSMFTReconstruction/Cluster.h"
#include <string>

using namespace o2::Base;
using o2::ITSMFT::Cluster;

void FromCluster2String(const Cluster &cluster, std::string &str);

#endif
