#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

typedef char           int8;
typedef unsigned char  uint8;
typedef short          int16;
typedef unsigned short uint16;
typedef int            int32;
typedef unsigned int   uint32;
typedef long           int64;
typedef unsigned long  uint64;

#define MAX_THREADS_PER_BLOCK 512

__global__ void
DijkstraKernel1(int* g_graph_nodes,
		int* g_graph_edges,
		int* g_graph_weights,
		int* g_up_cost,
		bool* g_graph_mask,
		int* g_cost,
		int no_of_nodes,
		int no_of_edges) {
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
DijkstraKernel2(int* g_graph_nodes,
		int* g_graph_edges,
		int* g_graph_weights, 
		int* g_up_cost,
		bool* g_graph_mask, 
		int* g_cost,
		bool *d_finished,
		int no_of_nodes,
		int no_of_edges) {
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

void copyToDevice(void *h_mem, void **d_mem, size_t size) {
  CUDA_CHECK_RETURN( cudaMalloc( d_mem, size) );
  CUDA_CHECK_RETURN( cudaMemcpy( *d_mem, h_mem, size, cudaMemcpyHostToDevice) );
}

int main( int argc, char** argv) {

  const int32 MAX_COST = 1 << 30;

  uint32 no_of_nodes = 0;
  uint32 no_of_edges = 0;
  
  printf("Reading File\n");
  // Read in Graph from a file
  FILE *fp = fopen("NY.out","r");
  if(!fp) {
    printf("Error Reading graph file\n");
    return -1;
  }

  int32 source = 0;

  fscanf(fp,"%d",&no_of_nodes);
  printf("No of Nodes: %d\n",no_of_nodes);

  // Make execution Parameters according to the number of nodes
  // Distribute threads across multiple Blocks if necessary
  uint32 num_of_blocks = 1;
  uint32 num_of_threads_per_block = no_of_nodes;
  if(no_of_nodes>MAX_THREADS_PER_BLOCK) {
    num_of_blocks = (uint32)ceil(no_of_nodes/(double)MAX_THREADS_PER_BLOCK);
    num_of_threads_per_block = MAX_THREADS_PER_BLOCK;
  }

  // allocate host memory
  int* h_graph_nodes = (int*) malloc(sizeof(int)*no_of_nodes);
  bool *h_graph_mask = (bool*) malloc(sizeof(bool)*no_of_nodes);
  int *h_up_cost = (int*) malloc(sizeof(int)*no_of_nodes);

  int start, edgeno;
  // initalize the memory
  for( uint32 i = 0; i < no_of_nodes; i++ ) {
    fscanf(fp,"%d %d",&start,&edgeno);
    h_graph_nodes[i] = start;
    h_up_cost[i] = MAX_COST;
    h_graph_mask[i]=false;
  }

  //read the source int from the file
  fscanf(fp,"%d",&source);
  printf("Source vid: %d\n", source);

  //set the source int as true in the mask
  h_graph_mask[source]=true;

  fscanf(fp,"%d",&no_of_edges);
  printf("No of Edges: %d\n", no_of_edges);

  int id;
  int* h_graph_edges = (int*) malloc(sizeof(int)*no_of_edges);
  int* h_graph_weights = (int*) malloc(sizeof(int)*no_of_edges);
  for(int i=0; i < no_of_edges ; i++) {
    fscanf(fp,"%d",&id);
    h_graph_edges[i] = id;
    fscanf(fp,"%d",&id);
    h_graph_weights[i] = id;
  }

  if(fp) fclose(fp);

  printf("Read File\n");
  printf("Avg Branching Factor: %f\n",no_of_edges/(float)no_of_nodes);

  // allocate mem for the result on host side
  int* h_cost = (int*) malloc( sizeof(int)*no_of_nodes);
  for(int i=0;i<no_of_nodes;i++) h_cost[i]= MAX_COST;
  h_cost[source]=0;

  int32 *d_cost, *d_up_cost, *d_graph_nodes, *d_graph_edges, *d_graph_weights;
  bool *d_graph_mask;

  // Copy lists to device memory
  copyToDevice( h_up_cost,       (void**)&d_up_cost,       sizeof(int32)*no_of_nodes );
  copyToDevice( h_cost,          (void**)&d_cost,          sizeof(int32)*no_of_nodes );
  copyToDevice( h_graph_nodes,   (void**)&d_graph_nodes,   sizeof(int32)*no_of_nodes );
  copyToDevice( h_graph_edges,   (void**)&d_graph_edges,   sizeof(int32)*no_of_edges );
  copyToDevice( h_graph_weights, (void**)&d_graph_weights, sizeof(int32)*no_of_edges );
  copyToDevice( h_graph_mask,	 (void**)&d_graph_mask,    sizeof(bool)*no_of_nodes );

  //make a bool to check if the execution is over
  bool finished, *d_finished;
  CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_finished, sizeof(bool)));

  // setup execution parameters
  dim3  grid( num_of_blocks, 1, 1);
  dim3  threads( num_of_threads_per_block, 1, 1);

  int k=0;
  do {
    DijkstraKernel1<<< grid, threads, 0 >>>( d_graph_nodes, d_graph_edges, d_graph_weights, d_up_cost,
					     d_graph_mask, d_cost, no_of_nodes, no_of_edges );
    k++;
    finished = false;
    CUDA_CHECK_RETURN( cudaMemcpy( d_finished, &finished, sizeof(bool), cudaMemcpyHostToDevice ) );
    DijkstraKernel2<<< grid, threads, 0 >>>( d_graph_nodes, d_graph_edges, d_graph_weights, d_up_cost,
					     d_graph_mask, d_cost, d_finished, no_of_nodes, no_of_edges );
    CUDA_CHECK_RETURN( cudaThreadSynchronize() );    // Wait for the GPU launched work to complete
    CUDA_CHECK_RETURN( cudaGetLastError() );
    CUDA_CHECK_RETURN( cudaMemcpy( &finished, d_finished, sizeof(bool), cudaMemcpyDeviceToHost ) );
  } while( finished );

  // copy result from device to host
  CUDA_CHECK_RETURN( cudaMemcpy( h_cost, d_cost, sizeof(int)*no_of_nodes, cudaMemcpyDeviceToHost) );

  printf("Computation finished\n");

  //Store the result into a file
  FILE *fpo = fopen("result.txt","w");
  for(int i=0;i<no_of_nodes;i++)
    fprintf(fpo,"%d) cost:%d\n",i,h_cost[i]);
  fclose(fpo);
  printf("Result stored in result.txt\n");

  // cleanup memory
  free(h_graph_nodes);
  free(h_graph_edges);
  free(h_graph_mask);
  free(h_graph_weights);
  free(h_up_cost);
  free(h_cost);
  CUDA_CHECK_RETURN(cudaFree(d_graph_nodes));
  CUDA_CHECK_RETURN(cudaFree(d_graph_edges));
  CUDA_CHECK_RETURN(cudaFree(d_graph_mask));
  CUDA_CHECK_RETURN(cudaFree(d_graph_weights));
  CUDA_CHECK_RETURN(cudaFree(d_up_cost));
  CUDA_CHECK_RETURN(cudaFree(d_cost));
  CUDA_CHECK_RETURN(cudaFree(d_finished));
  CUDA_CHECK_RETURN(cudaDeviceReset());
  return 0;
}
