#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "types.h"
#include "qt.hpp"

using namespace std;

#define MAX_THREADS_PER_BLOCK 512

__global__ void
DijkstraKernel1(uint32  no_of_nodes,
                uint32  no_of_edges,
                uint32 *g_graph_nodes,
                uint32 *g_graph_edges,
                uint32 *g_graph_weights,
                uint32 *g_up_cost,
                char   *g_graph_mask,
                uint32 *g_cost) {
  int tid = blockIdx.x*MAX_THREADS_PER_BLOCK + threadIdx.x;
  int i,end,id;
  if(tid<no_of_nodes && g_graph_mask[tid]) {
    end = (tid < no_of_nodes-1) ? g_graph_nodes[tid+1] : no_of_edges;
    for(i = g_graph_nodes[tid]; i< end; i++) {
      id = g_graph_edges[i];
      atomicMin(&g_up_cost[id], g_cost[tid]+g_graph_weights[i]);
    }
    g_graph_mask[tid]=false;
  }
}

__global__ void
DijkstraKernel2(uint32  no_of_nodes,
                uint32 *g_up_cost,
                char   *g_graph_mask,
                uint32 *g_cost,
                char   *d_finished) {
  int tid = blockIdx.x*MAX_THREADS_PER_BLOCK + threadIdx.x;
  if(tid<no_of_nodes && g_cost[tid] > g_up_cost[tid]) {
    g_cost[tid] = g_up_cost[tid];
    g_graph_mask[tid] = true;
    *d_finished = true;
  }
  if(tid<no_of_nodes) g_up_cost[tid] = g_cost[tid];
}

/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#define CUDA_CHECK_RETURN(value) {                                    \
    cudaError_t _m_cudaStat = value;                                  \
    if (_m_cudaStat != cudaSuccess) {                                 \
        fprintf(stderr, "Error %s at line %d in file %s\n",           \
                cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); \
        exit(1);                                                      \
    } }

struct device_ptrs {
  uint32 *up_cost;
  uint32 *cost;
  uint32 *nodes;
  uint32 *edges;
  uint32 *evals;
  char   *mask;
  char   *finished;
};

struct host_ptrs {
  uint32  nn;
  uint32  ne;
  uint32 *graph_nodes;
  uint32 *up_cost;
  uint32 *graph_edges;
  uint32 *graph_weights;
  char   *mask;
};

void graph_to_dev(device_ptrs dps,
                  uint32 nn, uint32 *nodes,
                  uint32 ne, uint32 *edges, uint32 *evals) {
  CUDA_CHECK_RETURN( cudaMemcpy( dps.nodes, nodes, 4*nn, cudaMemcpyHostToDevice) );
  CUDA_CHECK_RETURN( cudaMemcpy( dps.edges, edges, 4*ne, cudaMemcpyHostToDevice) );
  CUDA_CHECK_RETURN( cudaMemcpy( dps.evals, evals, 4*ne, cudaMemcpyHostToDevice) );
}

void prob_to_dev(device_ptrs dps,
                 uint32 nn, uint32 *up_cost, uint32 *cost, char *mask) {
  CUDA_CHECK_RETURN( cudaMemcpy( dps.up_cost, up_cost, 4*nn, cudaMemcpyHostToDevice) );
  CUDA_CHECK_RETURN( cudaMemcpy( dps.cost,    cost,    4*nn, cudaMemcpyHostToDevice) );
  CUDA_CHECK_RETURN( cudaMemcpy( dps.mask,    mask,    nn, cudaMemcpyHostToDevice) );
}

void sssp( device_ptrs dps, host_ptrs hps, const uint32 source_id,
           uint32 *h_cost) {
  cout << " sssp(" << source_id << ")..."; cout.flush();;
  const uint32 MAX_COST = 1 << 30;

  for( uint32 i=0; i<hps.nn; i++) {
    hps.up_cost[i] = MAX_COST;
    h_cost[i]      = MAX_COST;
    hps.mask[i]    = false;
  }

  h_cost[source_id] = 0;
  hps.mask[source_id] = true;

  // Copy lists to device memory
  prob_to_dev(dps, hps.nn, hps.up_cost, h_cost, hps.mask);

  //make a char to check if the execution is over
  char finished;

  // setup execution parameters
  // Make execution Parameters according to the number of nodes
  // Distribute threads across multiple Blocks if necessary
  uint32 num_of_blocks = 1;
  uint32 num_of_threads_per_block = hps.nn;
  if(hps.nn>MAX_THREADS_PER_BLOCK) {
    num_of_blocks = (uint32)ceil(hps.nn/(double)MAX_THREADS_PER_BLOCK);
    num_of_threads_per_block = MAX_THREADS_PER_BLOCK;
  }
  dim3  grid( num_of_blocks, 1, 1);
  dim3  threads( num_of_threads_per_block, 1, 1);

  cout << "running kernels..."; cout.flush();
  uint32 k=0;
  do {
    DijkstraKernel1<<< grid, threads, 0 >>>( hps.nn,
                                             hps.ne,
                                             dps.nodes,
                                             dps.edges,
                                             dps.evals,
                                             dps.up_cost,
                                             dps.mask,
                                             dps.cost );
    k++;
    CUDA_CHECK_RETURN( cudaMemset( dps.finished, 0, 1 ) );
    DijkstraKernel2<<< grid, threads, 0 >>>( hps.nn,
                                             dps.up_cost,
                                             dps.mask,
                                             dps.cost,
                                             dps.finished);
    CUDA_CHECK_RETURN( cudaThreadSynchronize() );    // Wait for the GPU launched work to complete
    CUDA_CHECK_RETURN( cudaGetLastError() );
    CUDA_CHECK_RETURN( cudaMemcpy( &finished, dps.finished, 1, cudaMemcpyDeviceToHost ) );
  } while( finished );
  cout << "done in " << k << " iterations..."; cout.flush();
  // copy result from device to host
  CUDA_CHECK_RETURN( cudaMemcpy( h_cost, dps.cost, 4*hps.nn, cudaMemcpyDeviceToHost) );
  cout << "done" << endl;
}

void decompose(const Qt &qt, const qblck b, vector<qblck> &l) {
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

void enque_allpairs(const Qt &qt, const vector<qblck> &as, const vector<qblck> &bs, workq &Q) {
  for(vector<qblck>::const_iterator a = as.begin(); a != as.end(); a++) {
    for(vector<qblck>::const_iterator b = bs.begin(); b != bs.end(); b++) {
      if(((*a) != (*b)) || (qt.isnotleaf(*a))) {
        Q.push_back(make_pair(*a, *b));
      }
    }
  }
}

struct qblck_info {
  typedef typename std::map<qblck, uint32> distmap;
  uint32   netdiam;
  Qvtx    *p;
  distmap *dists;
  
  qblck_info(uint32 d, Qvtx *pt) : netdiam(d), p(pt) { dists = new distmap; };
  ~qblck_info() { cout << "~" << endl; delete dists; };
};

void process_allpairs(const double sep, ofstream &rf,
                      const Qt &qt, vector<qblck> &as, vector<qblck> &bs,
                      workq &Q, device_ptrs dps, host_ptrs hps, uint32 *h_cost) {
  map<qblck, qblck_info*> qblock_cache;
  for(vector<qblck>::iterator a = as.begin(); a != as.end(); a++) {
    Qvtx *pa = qt.getRep(*a);
    sssp(dps, hps, pa->vid, h_cost);
    uint32 d = qt.netdiam(h_cost, *a);
    qblock_cache.insert(make_pair(*a, new qblck_info(d, pa)));

    for(vector<qblck>::iterator b = bs.begin(); b != bs.end(); b++) {
      if(((*a) == (*b)) && (qt.isnotleaf(*a))) { // Maintain vectors for indexing
      } else {
        // a and b distinct, lets compute sssps all around
        Qvtx *pb = qt.getRep(*b);
	qblock_cache.find(*a)->second->dists->insert(make_pair(*b, h_cost[pb->vid]));
      }
    }
  }
  for(vector<qblck>::iterator a = as.begin(); a != as.end(); a++) {
    qblck_info *qblcki = qblock_cache.find(*a)->second;
    for(uint64 i=0; i<bs.size(); i++) {
      const qblck *b = &bs.at(i);
      if(((*a) == (*b)) && (qt.isnotleaf(*a))) {
        Q.push_back(make_pair(*a, *b));
      } else {
        // a and b distinct, lets compute sssps all around
        Qvtx *pb      = qt.getRep(*b);
        uint32 dg_a_b = qblcki->dists->find(*b)->second;
	map<qblck, qblck_info*>::iterator bi = qblock_cache.find(*b);
	uint32 db = 0;
	if(bi != qblock_cache.end()) {
	  db = bi->second->netdiam;
	} else {
	  sssp(dps, hps, pb->vid, h_cost);
	  // Measure diameter of B
	  db = qt.netdiam(h_cost, *b);
	}
        // r = max(da, db)
        uint32 r = std::max(qblcki->netdiam, db);
        // if dg/r >= sep
        if( dg_a_b >= sep*r ) {
          cout << " L: " << hex << CODE_OF_QBLCK(*a) << " -> ";
          cout << CODE_OF_QBLCK(*b) << " = " << dec << dg_a_b << endl;
          rf << CODE_OF_QBLCK(*a) << " -> " << CODE_OF_QBLCK(*b) << " = " << dg_a_b << endl;
        } else {
          cout << " ~L" << endl;
          std::vector<qblck> la, lb;
          decompose(qt, *a, la);
          decompose(qt, *b, lb);
          enque_allpairs(qt, la, lb, Q);
        }
      }
    }
  }
}

int build_oracle() {

  uint32 nn = 0;
  uint32 ne = 0;

  printf("Reading File\n");
  // Read in Graph from a file
  FILE *fp = fopen("/home/natep/cuda-workspace/sssp/NY.out","r");
  if(!fp) {
    printf("Error Reading graph file\n");
    return -1;
  }

  int32 source_id = 0;

  fscanf(fp,"%d",&nn);
  printf("No of Nodes: %d\n",nn);

  // allocate host memory
  uint32 *h_graph_nodes, *h_up_cost;
  char   *h_mask;
  CUDA_CHECK_RETURN( cudaMallocHost( &h_graph_nodes, sizeof(uint32)*nn ) );
  CUDA_CHECK_RETURN( cudaMallocHost( &h_up_cost,     sizeof(uint32)*nn ) );
  CUDA_CHECK_RETURN( cudaMallocHost( &h_mask,        nn ) );

  // initalize the memory
  uint32 start, edgeno;
  for( uint32 i = 0; i < nn; i++ ) {
    fscanf(fp,"%d %d",&start,&edgeno);
    h_graph_nodes[i] = start;
  }

  //read the source int from the file
  fscanf(fp,"%d",&source_id);
  printf("Source vid: %d\n", source_id);

  fscanf(fp,"%d",&ne);
  printf("No of Edges: %d\n", ne);

  uint32 id;
  uint32* h_graph_edges, *h_graph_weights;
  CUDA_CHECK_RETURN( cudaMallocHost( &h_graph_edges, sizeof(uint32)*ne ) );
  CUDA_CHECK_RETURN( cudaMallocHost( &h_graph_weights, sizeof(uint32)*ne ) );
  for(uint32 i=0; i < ne ; i++) {
    fscanf(fp,"%d",&id);
    h_graph_edges[i] = id;
    fscanf(fp,"%d",&id);
    h_graph_weights[i] = id;
  }

  if(fp) fclose(fp);

  // Put all host pointers (except h_cost) together
  host_ptrs hps;
  hps.nn = nn;
  hps.ne = ne;
  hps.graph_nodes = h_graph_nodes;
  hps.up_cost = h_up_cost;
  hps.graph_edges = h_graph_edges;
  hps.graph_weights = h_graph_weights;
  hps.mask = h_mask;

  printf("Read File\n");
  printf("Avg Branching Factor: %f\n",ne/(float)nn);

  // allocate mem for the result on host side
  uint32 *h_cost;
  CUDA_CHECK_RETURN( cudaMallocHost( &h_cost, sizeof(uint32)*nn ) );

  // allocate everything in device
  device_ptrs dps;
  CUDA_CHECK_RETURN( cudaMalloc( &dps.up_cost, 4*nn ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.cost,    4*nn ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.nodes,   4*nn ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.edges,   4*ne ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.evals,   4*ne ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.mask,    nn ) );
  CUDA_CHECK_RETURN( cudaMalloc( &dps.finished,1));
  graph_to_dev( dps,
                nn, h_graph_nodes,
                ne, h_graph_edges, h_graph_weights);

  /********************************************************************************/

  ifstream cof("/home/natep/cuda-workspace/sssp/NY.co");
  string v;
  uint32 lat, lon;
  uint32 max = 1U << 28;
  vector<Qvtx*> qvtxes;
  qvtxes.reserve(nn);
  Qt qt;
  for(uint32 i=0; i<nn; i++) {
    cof >> v >> id >> lat >> lon;
    qvtxes.push_back(new Qvtx(i, morton_code(max+lat, max+lon)));
    qt.insert(qvtxes.at(i));
  }
  cof.close();

  printf("all read, first=%lu, qt.size=%lu\n", qvtxes[0]->z, qt.size());

  /********************************************************************************/

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
        process_allpairs(sep, rf, qt, l, l, Q, dps, hps, h_cost);
      } else { // do nothing
        cout << " Same leaf" << endl;
      }
    } else {
      cout << " a!=b" << endl;
      // Choose rep point of A
      Qvtx *pa = qt.getRep(a);
      if(NULL == pa) {
        cout << "nonexistent node ended up in q as a" << endl;
        continue;
      }
      // Get sssp from pa
      sssp(dps, hps, pa->vid, h_cost);
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
      sssp(dps, hps, pb->vid, h_cost);
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
        process_allpairs(sep, rf, qt, la, lb, Q, dps, hps, h_cost);
      }
    }
  }

  /********************************************************************************/

  printf("Computation finished\n");

  //Store the result into a file
  rf.close();
  printf("Result stored in result.txt\n");

  // cleanup memory
  CUDA_CHECK_RETURN(cudaFreeHost(h_graph_nodes));
  CUDA_CHECK_RETURN(cudaFreeHost(h_graph_edges));
  CUDA_CHECK_RETURN(cudaFreeHost(h_graph_weights));
  CUDA_CHECK_RETURN(cudaFreeHost(h_mask));
  CUDA_CHECK_RETURN(cudaFreeHost(h_up_cost));
  CUDA_CHECK_RETURN(cudaFreeHost(h_cost));
  CUDA_CHECK_RETURN(cudaFree(dps.nodes));
  CUDA_CHECK_RETURN(cudaFree(dps.edges));
  CUDA_CHECK_RETURN(cudaFree(dps.mask));
  CUDA_CHECK_RETURN(cudaFree(dps.evals));
  CUDA_CHECK_RETURN(cudaFree(dps.up_cost));
  CUDA_CHECK_RETURN(cudaFree(dps.cost));
  CUDA_CHECK_RETURN(cudaFree(dps.finished));
  CUDA_CHECK_RETURN(cudaDeviceReset());
  return 0;
}

int32 main( int32 argc, char** argv) {
  return build_oracle();
}
