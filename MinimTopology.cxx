#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "./MinimTopology.h"

MinimTopology::MinimTopology():fPattern(),fHash(0){
}

MinimTopology::MinimTopology(const std::string &str) : fHash(0) {
  SetPattern(str);
}

void MinimTopology::SetPattern(const std::string &str) {
  int nBytes = (int)str.size();
  fPattern.resize(nBytes,0);
  memcpy(&fPattern[0],&str[0],nBytes);
  nBytes-=2;
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
    std::cout << "ERROR: no fired pixels\n";
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

std::ostream& operator<<(std::ostream& os, const MinimTopology& top){
	int rs = top.fPattern[0];
	int cs = top.fPattern[1];
	os << "rs: " << rs << " cs: " << cs << " #bytes: " << top.fPattern.length() << std::endl;
	unsigned char tempChar = 0;
	int s=0;
	int ic = 0;
  for (unsigned int i=2; i<top.fPattern.length(); i++){
		tempChar = top.fPattern[i];
		s=128; //0b10000000
    while(s>0){
			if(ic%cs==0) os << "|";
			ic++;
      if((tempChar&s)!=0) os << '+';
      else os << ' ';
      s/=2;
			if(ic%cs==0) os << "|" << std::endl;
			if(ic==(rs*cs)) break;
    }
		if(ic==(rs*cs)) break;
  }
  os<< std::endl;
  return os;
}
