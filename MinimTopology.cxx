#include <iostream>
#include <string>
#include "./MinimTopology.h"

using namespace std;

MinimTopology::MinimTopology():fPattern(),fHash(0){
}

MinimTopology::~MinimTopology(){
}

MinimTopology::MinimTopology(const AliITSMFTClusterPix &cluster) : fHash(0) {
  int rs = cluster.GetPatternRowSpan();
  int cs = cluster.GetPatternColSpan();
  fPattern.push_back(rs);
	fPattern.push_back(cs);
  unsigned char tempChar = 0;
  int BitCounter=7;
  for(int ir=0; ir<rs; ir++){
    for(int ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	      fPattern.push_back(tempChar);
	      tempChar=0;
	      BitCounter=7;
      }
      if(cluster.TestPixel(ir,ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  fPattern.push_back(tempChar);
  int nBytes = fPattern.length()-2;
  fHash = ((unsigned long)(hashFunction(fPattern.data(),fPattern.length())))<<32;
  if(nBytes>=4){
    fHash += ((((unsigned long)fPattern[2])<<24) + (((unsigned long)fPattern[3])<<16) + (((unsigned long)fPattern[4])<<8) + ((unsigned long)fPattern[5]));
  }
  else if(nBytes==3){
    fHash += ((((unsigned long)fPattern[2])<<24) + (((unsigned long)fPattern[3])<<16) + (((unsigned long)fPattern[4])<<8));
  }
  else if(nBytes==2){
    fHash += ((((unsigned long)fPattern[2])<<24) + (((unsigned long)fPattern[3])<<16));
  }
  else if(nBytes==1){
    fHash += ((((unsigned long)fPattern[2])<<24));
  }
  else{
    cout << "ERROR: no fired pixels\n";
    exit(1);
  }
}

unsigned int MinimTopology::hashFunction(const void* key, int len){
  //
  //Developed from https://github.com/rurban/smhasher , function MurMur2
  //
  // 'm' and 'r' are mixing constants generated offline.
  const unsigned int m =0x5bd1e995;
  const int r = 24;
  // Initialize the hash
  unsigned int h = len^0xdeadbeef;
  // Mix 4 bytes at a time into the hash
  const unsigned char* data = (const unsigned char *)key;
  //int recIndex=0;
  while(len >= 4){
    unsigned int k = *(unsigned int*)data;
    k *= m;
    k ^= k >> r;
    k *= m;
    h *= m;
    h ^= k;
    data += 4;
    len -= 4;
  }
  // Handle the last few bytes of the input array
  switch(len){
	  case 3: h ^= data[2] << 16;
	  case 2: h ^= data[1] << 8;
	  case 1: h ^= data[0];
	  h *= m;
  };
  // Do a few final mixes of the hash to ensure the last few
  // bytes are well-incorporated.
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  return h;
}

std::ostream& MinimTopology::printTop(std::ostream &out){
	int rs = fPattern[0];
	int cs = fPattern[1];
	out << "rs: " << rs << " cs: " << cs << " #bytes: " << fPattern.length() << endl;
	unsigned char tempChar = 0;
	int s=0;
	int ic = 0;
  for (unsigned int i=2; i<fPattern.length(); i++){
		tempChar = fPattern[i];
		s=128; //0b10000000
    while(s>0){
			if(ic%cs==0) out << "|";
			ic++;
      if((tempChar&s)!=0) out << '+';
      else out << ' ';
      s/=2;
			if(ic%cs==0) out << "|" << endl;
			if(ic==(rs*cs)) break;
    }
		if(ic==(rs*cs)) break;
  }
  out<< endl;
}
