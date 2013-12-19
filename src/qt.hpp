#ifndef _QT_H_
#define _QT_H_

#include <string.h>
#include <iostream>
#include <algorithm> // std::min
#include "types.h"

const uint64 MAXUINT64 = ~0;

typedef uint64 zcode;
typedef uint64 level;
typedef uint64 qblck;

// Create bitmask of w 1s
#define BITMASK(w) ((1UL << (w)) - 1UL)

// Create bitmask of w 0s
#define INVMASK(w) (~(BITMASK(w)))

#define HALF_WDTH 29
#define CODE_SIZE 58
#define CODE_MASK (BITMASK(CODE_SIZE))

#define CODE_OF_QBLCK(b) (b & CODE_MASK)
#define LEVEL_OF_QBLCK(b) (b >> CODE_SIZE)

#define QBLCK(l, c) ((static_cast<level>(l) << CODE_SIZE) | \
		     (c & (INVMASK(CODE_SIZE-2*(l)))))

#define COUT_QBLCK(b) \
    "qblck(" << b << "=" << LEVEL_OF_QBLCK(b) << "|" << CODE_OF_QBLCK(b) << ")"

inline zcode expand(uint64 x) {
  x = ( x | (x << 16)) & 0x0000ffff0000ffff;
  x = ( x | (x <<  8)) & 0x00ff00ff00ff00ff;
  x = ( x | (x <<  4)) & 0x0f0f0f0f0f0f0f0f;
  x = ( x | (x <<  2)) & 0x3333333333333333;
  x = ( x | (x <<  1)) & 0x5555555555555555;
  return x;
}

inline zcode morton_code(uint64 x, uint64 y) {
  // Interleave bits of x and y, so that the lower CODE_SIZE/2
  // bits of x are in the even positions and y in the odd;

  return expand(x) | (expand(y) << 1);
}

inline uint64 contract(uint64 x) {
  x = ( x            ) & 0x5555555555555555;
  x = ( x | (x >>  1)) & 0x3333333333333333;
  x = ( x | (x >>  2)) & 0x0f0f0f0f0f0f0f0f;
  x = ( x | (x >>  4)) & 0x00ff00ff00ff00ff;
  x = ( x | (x >>  8)) & 0x0000ffff0000ffff;
  x = ( x | (x >> 16)) & 0x00000000ffffffff;
  return x;
}

inline latlon morton_uncode(zcode z) {
  uint64 x = contract(z);
  uint64 y = contract(z>>1);
  return std::make_pair(x, y);
}

inline qblck child(qblck b, uint64 c) {
  level l = LEVEL_OF_QBLCK(b)+1;
  qblck bc = QBLCK(l, ((CODE_OF_QBLCK(b)) | (c << (CODE_SIZE-2*l))));
  return bc;
}

inline qblck child00(qblck b) { return child(b, 0); }

inline qblck child01(qblck b) { return child(b, 1); }

inline qblck child10(qblck b) { return child(b, 2); }

inline qblck child11(qblck b) { return child(b, 3); }

struct Qvtx {
  uint64 vid;
  zcode z;
  Qvtx(uint64 vid, zcode z)
    : vid(vid), z(z) {}
  Qvtx(uint64 vid, uint32 x, uint32 y)
    : vid(vid), z(morton_code(x, y)) {}
  Qvtx(const Qvtx &v)
    : vid(v.vid), z(v.z) {}
  ~Qvtx() {}
};

//  qblck, Qvtx
struct qtree_iterator {
  void increment();
  qblck qblock();
  Qvtx* qvtx();
};

struct Qt {
  Qt();
  ~Qt();
  bool contains(qblck b);
  void insert(Qvtx *v);
  uint64 size();
  Qvtx* getRep(qblck b);
  qtree_iterator block_it(qblck b);
};

// -1 if l <z r, 1 if r <z l, 0 if one contains the other
int cmp_qblck(const qblck &l, const qblck &r) {
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

#endif // _QT_H_
