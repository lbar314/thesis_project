#ifndef CLUSTER2STRING_H
#define CLUSTER2STRING_H

#include "AliITSMFTClusterPix.h"
#include <string>

void FromCluster2String(const AliITSMFTClusterPix &cluster, std::string &str);

#endif
