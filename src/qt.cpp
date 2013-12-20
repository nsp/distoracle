#include "qt.hpp"
#include "btree_map.h"


struct QblckComparer
  : public btree::btree_key_compare_to_tag {
  int operator()(const qblck &l, const qblck &r) const {
    return cmp_qblck(l, r);
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

bool Qt::contains(const qblck b) const {
  Qmap::iterator lookup = qmap->find(b);
  return lookup != qmap->end();
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

uint64 Qt::size() const {
  return qmap->size();
}

bool Qt::isleaf(const qblck b) const {
  Qmap::const_iterator lookup = qmap->find(b);
  if(lookup == qmap->end()) {
    return false;
  }
  return lookup.key() == b;
}

bool Qt::isnotleaf(const qblck b) const {
  return !isleaf(b);
}

uint64 crowdist2(const latlon a, const latlon b) {
  uint64 xd = std::max(a.first, b.first) - std::min(a.first, b.first);
  uint64 yd = std::max(a.second, b.second) - std::min(a.second, b.second);
  return xd*xd + yd*yd;
}

Qvtx* Qt::getRep(const qblck b) const {
  Qmap::const_iterator lookup = qmap->find(b);
  if(lookup == qmap->end()) { // no block
    return NULL;
  } else if(lookup.key() == b) { // is leaf
    return lookup->second;
  }
  // Nonleaf, return closest
  latlon cxy = morton_uncode(CODE_OF_QBLCK(child11(b)));
  Qvtx *minvtx = lookup->second;
  uint64 mindst = crowdist2(morton_uncode(minvtx->z), cxy);
  uint64 dst;
  while(++lookup != qmap->end() && cmp_qblck(lookup.key(), b) == 0) {
    zcode oz = lookup->second->z;
    latlon oxy = morton_uncode(oz);
    dst = crowdist2(oxy, cxy);
    if(dst < mindst) {
      minvtx = lookup->second;
      mindst = dst;
    }
  }
  return minvtx;
}

uint32 Qt::netdiam(const uint32 *const dists, const qblck b) const {
  if(isleaf(b)) {
    return 0;
  }
  uint32 maxdst = 0, dst;
  Qmap::const_iterator lookup = qmap->find(b);
  while(lookup != qmap->end() && cmp_qblck(lookup.key(), b) == 0) {
    dst = dists[lookup->second->vid];
    if(dst >= maxdst) {
      maxdst = dst;
    }
    lookup++;
  }
  return maxdst;
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

using namespace std;

void decompose(Qt &qt, const qblck b, vector<qblck> &l) {
  if(qt.isleaf(b)) {
    l.push_back(b);
  } else {
    for(uint64 d=0; d<4; d++) {
      qblck c = child(b, d);
      if(qt.contains(c)) {
	l.push_back(c);
      }
    }
  }
}

void process_allpairs(const Qt &qt, const vector<qblck> &as, const vector<qblck> &bs, workq &Q) {
  for(vector<qblck>::const_iterator a = as.begin(); a != as.end(); a++) {
    for(vector<qblck>::const_iterator b = bs.begin(); b != bs.end(); b++) {
      if(((*a) != (*b)) || (qt.isnotleaf(*a))) {
	Q.push_back(make_pair(*a, *b));
      }
    }
  }
}

int main(int argc, char* argv[]) {
  
  uint64 nn = 264346;
  std::ifstream cof("/home/natep/cuda-workspace/sssp/NY.co");
  std::string v;
  int32 id, lat, lon;
  std::vector<Qvtx*> qvtxes;
  qvtxes.reserve(nn);
  Qt qt;
  for(uint32 i=0; i<4; i++) {
    cof >> v >> id >> lat >> lon;
    qvtxes.push_back(new Qvtx(i, morton_code(std::abs(lat), lon)));
    qt.insert(qvtxes.at(i));
  }
  cof.close();

  printf("all read, qt.size=%lu\n", qt.size());

  uint32 *h_cost = (uint32*)malloc(nn*sizeof(uint32));
  //Store the result into a file
  ofstream rf("/home/natep/cuda-workspace/sssp/result.txt");
  double eps = 0.5;
  double sep = 2/eps;
  std::deque<std::pair<qblck, qblck> > Q;
  qblck root = QBLCK(0, 0);
  Q.push_back(make_pair(root, root));
  uint64 qiters = 0;
  while(!Q.empty()) {
    qblck a = Q.front().first;
    qblck b = Q.front().second;
    Q.pop_front();
    cout << ++qiters << "/" << Q.size() << ": ";
    cout << "a=" << LEVEL_OF_QBLCK(a) << "|" << hex << CODE_OF_QBLCK(a) << dec << ", ";
    cout << "b=" << LEVEL_OF_QBLCK(b) << "|" << hex << CODE_OF_QBLCK(b) << dec << endl;
    cerr << qiters << endl;
    if(a==b) {
      if(qt.isnotleaf(a)) {
	cout << " Same nonleaf" << endl;
	vector<qblck> l;
	decompose(qt, a, l);
	process_allpairs(qt, l, l, Q);
      } else {
	cout << " Same leaf" << endl;
      }
    } else { // do nothing
      // Choose rep point of A
      Qvtx *pa = qt.getRep(a);
      if(NULL == pa) {
        cout << "nonexistent node ended up in q as a" << endl;
        continue;
      }
      // Get sssp from pa
      cout << " sssp(" << pa->vid << "|" << pa->z <<")..."; cout.flush();;
      // sssp(dps,
      //      nn, h_graph_nodes, h_up_cost, h_mask,
      //      ne, h_graph_edges, h_graph_weights,
      //      pa->vid,
      //      h_cost);
      cout << "done" << endl;
      // Measure diameter of A
      uint32 da = qt.netdiam(h_cost, a);
      // Choose rep point of B
      Qvtx *pb = qt.getRep(b);
      if(NULL == pa) {
        cout << "nonexistent node ended up in q as b" << endl;
        continue;
      }
      // dg = graph_dist(pa, pb)
      uint32 dg_a_b = h_cost[pb->vid];
      // Get sssp from pb
      cout << " sssp(" << pb->vid << "|" << pb->z <<")..."; cout.flush();;
      // sssp(dps,
      //      nn, h_graph_nodes, h_up_cost, h_mask,
      //      ne, h_graph_edges, h_graph_weights,
      //      pb->vid,
      //      h_cost);
      cout << "done" << endl;
      // Measure diameter of B
      uint32 db = qt.netdiam(h_cost, b);
      // r = max(da, db)
      uint32 r = std::max(da, db);
      // if dg/r >= sep
      if( dg_a_b >= sep*r ) {
        cout << " L: " << hex << CODE_OF_QBLCK(a) << " -> " << CODE_OF_QBLCK(b) << " = " << dec << dg_a_b << endl;
        rf << CODE_OF_QBLCK(a) << " -> " << CODE_OF_QBLCK(b) << " = " << dg_a_b << endl;
      } else {
        cout << " ~L" << endl;
        std::vector<qblck> la, lb;
	decompose(qt, a, la);
	decompose(qt, b, lb);
	process_allpairs(qt, la, lb, Q);
      }
    }
  }
  
  /********************************************************************************/

  std::cout << "Done" << std::endl;
}

#endif
