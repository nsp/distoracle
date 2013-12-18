#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_THREADS_PER_BLOCK 512
#define MAX_COST 10000000

int no_of_nodes=0;
int edge_list_size=0;
FILE *fp;

__global__ void
DijkastraKernel(int* g_graph_nodes, int* g_graph_edges,short int* g_graph_weights,
				int* g_graph_updating_cost, bool* g_graph_mask,
				int* g_cost , int no_of_nodes, int edge_list_size) {
	int tid = blockIdx.x*MAX_THREADS_PER_BLOCK + threadIdx.x;
	int i,end,id;
	if(tid<no_of_nodes && g_graph_mask[tid]) {
		if(tid < no_of_nodes-1)
			end = g_graph_nodes[tid+1];
		else
			end = edge_list_size;
		for(i = g_graph_nodes[tid]; i< end; i++) {
			id = g_graph_edges[i];
			atomicMin(&g_graph_updating_cost[id], g_cost[tid]+g_graph_weights[i]);
		}
		g_graph_mask[tid]=false;
	}
}


/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

int main( int argc, char** argv) {

    printf("Reading File\n");
    //Read in Graph from a file
    fp = fopen("samplegraph","r");
    if(!fp) {
      printf("Error Reading graph file\n");
      return -1;
    }

    int source = 0;

    fscanf(fp,"%d",&no_of_nodes);
    printf("No of Nodes: %d\n",no_of_nodes);

    int num_of_blocks = 1;
    int num_of_threads_per_block = no_of_nodes;

    //Make execution Parameters according to the number of nodes
    //Distribute threads across multiple Blocks if necessary
    if(no_of_nodes>MAX_THREADS_PER_BLOCK) {
      num_of_blocks = (int)ceil(no_of_nodes/(double)MAX_THREADS_PER_BLOCK);
      num_of_threads_per_block = MAX_THREADS_PER_BLOCK;
    }

    // allocate host memory
    int* h_graph_nodes = (int*) malloc(sizeof(int)*no_of_nodes);
    bool *h_graph_mask = (bool*) malloc(sizeof(bool)*no_of_nodes);
    int *h_graph_updating_cost = (int*) malloc(sizeof(int)*no_of_nodes);

    int start, edgeno;
    // initalize the memory
    int no=0;
    for( unsigned int i = 0; i < no_of_nodes; i++) {
      fscanf(fp,"%d %d",&start,&edgeno);
      if(edgeno>100)
	no++;
      h_graph_nodes[i] = start;
      h_graph_updating_cost[i] = MAX_COST;
      h_graph_mask[i]=false;
    }

    //read the source int from the file
    fscanf(fp,"%d",&source);

    //set the source int as true in the mask
    h_graph_mask[source]=true;
    //h_graph_counter[source]=0;

    fscanf(fp,"%d",&edge_list_size);

    int id;
    int* h_graph_edges = (int*) malloc(sizeof(int)*edge_list_size);
    short int* h_graph_weights = (short int*) malloc(sizeof(short int)*edge_list_size);
    for(int i=0; i < edge_list_size ; i++) {
      fscanf(fp,"%d",&id);
      h_graph_edges[i] = id;
      fscanf(fp,"%d",&id);
      h_graph_weights[i] = id;
    }

    if(fp) fclose(fp);

    printf("Read File\n");
    printf("Total %d dense nodes, Avg Branching Factor: %f\n",no,edge_list_size/(float)no_of_nodes);

    //Copy the int list to device memory
    int* d_graph_nodes;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_graph_nodes, sizeof(int)*no_of_nodes) );
    CUDA_CHECK_RETURN( cudaMemcpy( d_graph_nodes, h_graph_nodes, sizeof(int)*no_of_nodes, cudaMemcpyHostToDevice) );

    //Copy the Edge List to device Memory
    int* d_graph_edges;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_graph_edges, sizeof(int)*edge_list_size) );
    CUDA_CHECK_RETURN( cudaMemcpy( d_graph_edges, h_graph_edges, sizeof(int)*edge_list_size, cudaMemcpyHostToDevice) );

    short int* d_graph_weights;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_graph_weights, sizeof(short int)*edge_list_size) );
    CUDA_CHECK_RETURN( cudaMemcpy( d_graph_weights, h_graph_weights, sizeof(short int)*edge_list_size, cudaMemcpyHostToDevice) );

    //Copy the Mask to device memory
    bool* d_graph_mask;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_graph_mask, sizeof(bool)*no_of_nodes) );
    CUDA_CHECK_RETURN( cudaMemcpy( d_graph_mask, h_graph_mask, sizeof(bool)*no_of_nodes, cudaMemcpyHostToDevice) );

    // allocate mem for the result on host side
    int* h_cost = (int*) malloc( sizeof(int)*no_of_nodes);
    for(int i=0;i<no_of_nodes;i++) h_cost[i]= MAX_COST;
    h_cost[source]=0;
    // allocate device memory for result
    int* d_cost;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_cost, sizeof(int)*no_of_nodes));
    CUDA_CHECK_RETURN( cudaMemcpy( d_cost, h_cost, sizeof(int)*no_of_nodes, cudaMemcpyHostToDevice) );

    int* d_graph_updating_cost;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_graph_updating_cost, sizeof(int)*no_of_nodes));
    CUDA_CHECK_RETURN( cudaMemcpy( d_graph_updating_cost, h_graph_updating_cost, sizeof(int)*no_of_nodes, cudaMemcpyHostToDevice) );

    //make a bool to check if the execution is over

    bool *d_finished;
    bool finished;
    CUDA_CHECK_RETURN( cudaMalloc( (void**) &d_finished, sizeof(bool)));

    // setup execution parameters
    dim3  grid( num_of_blocks, 1, 1);
    dim3  threads( num_of_threads_per_block, 1, 1);

    DijkastraKernel<<< grid, threads, 0 >>>( d_graph_nodes, d_graph_edges, d_graph_weights, d_graph_updating_cost,
    		d_graph_mask, d_cost, no_of_nodes, edge_list_size);
    CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
    CUDA_CHECK_RETURN(cudaGetLastError());
    CUDA_CHECK_RETURN( cudaMemcpy( d_finished, &finished, sizeof(bool), cudaMemcpyHostToDevice) );

    // copy result from device to host
    CUDA_CHECK_RETURN( cudaMemcpy( h_cost, d_cost, sizeof(int)*no_of_nodes, cudaMemcpyDeviceToHost) );


    //Store the result into a file
    FILE *fpo = fopen("result.txt","w");
    for(int i=0;i<no_of_nodes;i++)
      fprintf(fpo,"%d) cost:%d\n",i,h_cost[i]);
    fclose(fpo);
    printf("Result stored in result.txt\n");

    // cleanup memory
    free( h_graph_nodes);
    free( h_graph_edges);
    free( h_graph_mask);
    free( h_graph_weights);
    free( h_graph_updating_cost);
    free( h_cost);
    CUDA_CHECK_RETURN(cudaFree(d_graph_nodes));
    CUDA_CHECK_RETURN(cudaFree(d_graph_edges));
    CUDA_CHECK_RETURN(cudaFree(d_graph_mask));
    CUDA_CHECK_RETURN(cudaFree(d_graph_weights));
    CUDA_CHECK_RETURN(cudaFree(d_graph_updating_cost));
    CUDA_CHECK_RETURN(cudaFree(d_cost));
    CUDA_CHECK_RETURN(cudaFree(d_finished));
	CUDA_CHECK_RETURN(cudaDeviceReset());
	return 0;
}
