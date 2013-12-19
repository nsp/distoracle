#include "qt.hpp"
#include "btree_map.h"

struct QblckComparer
  : public btree::btree_key_compare_to_tag {
  int operator()(const qblck &l, const qblck &r) const {
    level ll = LEVEL_OF_QBLCK(l);
    level rl = LEVEL_OF_QBLCK(r);
    zcode lz = CODE_OF_QBLCK(l);
    zcode rz = CODE_OF_QBLCK(r);
    // zcode diff = 0;
    if( ll == rl ) {
      // diff = lz - rz;
    } else if( ll < rl ) { // demote r and compare
      zcode oz = CODE_OF_QBLCK(QBLCK(ll, rz));
      // diff = lz - oz;
      rz = oz;
    } else { // demote l and compare
      zcode oz = CODE_OF_QBLCK(QBLCK(rl, lz));
      //diff = oz - rz;
      lz = oz;
    }
    // return (0 < diff) - (diff < 0);
    return (lz < rz) ? -1 : ((lz > rz) ? 1 : 0);
  }
};

typedef btree::btree_map<qblck, Qvtx*, QblckComparer> Qmap;

Qmap *qmap;

Qt::Qt() {
  qmap = new Qmap;
}

Qt::~Qt() {
  delete qmap;
}

bool Qt::contains(qblck b) {
  return false;
}

void Qt::insert(Qvtx *v) {
  qblck b = QBLCK(29, v->z);
  Qmap::iterator lookup = qmap->find(b);
  // constaining block already exists, retrieve and re-add
  if(lookup != qmap->end()) {
    qblck ob = lookup->first;
    level ol = LEVEL_OF_QBLCK(ob);
    Qvtx *ov = lookup->second;
    if(v->z == ov->z) { return; } // Can't handle identical points
    qmap->erase(lookup);
    for( level l=ol; l<=29; l++ ) {
      b = QBLCK(l, v->z);
      ob = QBLCK(l, ov->z);
      if(QblckComparer()(b, ob) != 0) {
	qmap->insert(std::make_pair(b, v));
	qmap->insert(std::make_pair(ob, ov));
	return;
      }
    }
    std::cerr << "max level" << std::endl;
    exit(22);
  } else { // no block exists, keep making smaller blocks
    for( level l=0; l<=29; l++ ) {
      b = QBLCK(l, v->z);
      lookup = qmap->find(b);
      if(lookup == qmap->end()) { // block does not exist
	  qmap->insert(std::make_pair(b, v));
	  return;
      }
    }
    std::cerr << "max level" << std::endl;
    exit(23); // shouldn't happen
  }
}

uint64 Qt::size() {
  return qmap->size();
}



Qvtx* Qt::getRep(qblck b) {
  pair<uint64, uint64> cxy = mo
  uint64 dir = rand() % 4;
  return NULL;
}

#if 0

#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
  // uint64 x = 0x48392;
  // uint64 y = 0x72390;
  // zcode z = morton_code(x, y);
  // std::pair<uint64, uint64> p = morton_uncode(z);
  // uint64 px = p.first;
  // uint64 py = p.second;
  // printf("%lu %lu %lu %lu %lu\n", x, y, z, px, py);
  //                    0123456701234567
  // qblck r = QBLCK(28, 0x03ffffffffffffff);
  // std::cout << std::hex << (1 << (CODE_SIZE-2-(2*28))) << std::endl;
  // qblck c00 = child00(r);
  // qblck c01 = child01(r);
  // qblck c10 = child10(r);
  // qblck c11 = child11(r);

  // std::cout << std::hex <<
  //   COUT_QBLCK(r) << " " << std::endl << 
  //   COUT_QBLCK(c00) << " " << std::endl << 
  //   COUT_QBLCK(c01) << " " << std::endl << 
  //   COUT_QBLCK(c10) << " " << std::endl << 
  //   COUT_QBLCK(c11) << " " << std::endl;
  
  uint64 nn = 264346;
  std::ifstream cof("NY.co");
  std::string v;
  int32 id, lat, lon;
  std::vector<Qvtx> qvtxes;
  qvtxes.reserve(nn);
  Qt qt;
  for(uint32 i=0; i<nn; i++) {
    cof >> v >> id >> lat >> lon;
    qvtxes.push_back(Qvtx(i, morton_code(std::abs(lat), lon)));
    qt.insert(&qvtxes.at(i));
  }
  cof.close();

  printf("all read, qt.size=%lu\n", qt.size());

}

#endif
