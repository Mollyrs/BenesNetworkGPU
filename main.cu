//put C:/Users/molly/Desktop/289Q/project/main.cu
//nvcc -std=c++11 main.cu

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <cooperative_groups.h>

// includes, project
#include <cuda.h>
#include <cuda_runtime.h>

using namespace cooperative_groups;
namespace cg = cooperative_groups;
//constucting 8x8 benes network
//four rows, 5 columns, 20 routers total

__host__
void makeLUT(int N, int* LUT){
	int M = N;
	int even = 0;
	int odd = 1;
	int LUTsize = N*(log2((double)N)*2 - 2);
	for (int i =0; i < LUTsize/2; i+=N){
		for (int j=0; j<N; j+=M){
			for (int k =0; k<M/2; k++){
				LUT[i+j+k] = even;
				even+=2;
			}
			for (int k =M/2; k<M; k++){
				LUT[i+j+k] = odd;
				odd+=2;
			}
		} even=0; odd=1; M = M/2;
	}
	for (int x=LUTsize-N, i=LUTsize/2; i<LUTsize;i+=N, x-=N){
		for(int j=0; j<N; j++){
			int newIndex = LUT[x+j-LUTsize/2];
			LUT[newIndex + i] = j;
		}
	}
	return;
}

int createMask(int n)
{
   int r = 0;
   for (int i=0; i<n; i++)
       r |= 1 << i;

   return r;
}


__global__
void benes(int N,  char* network, int* LUT, volatile int* valid, int mask){
	int idx = threadIdx.x;
	int in1, in2, in1_index, in2_index;
	// int level = blockIdx.x;

	// auto g = this_thread_block();
	// thread_group tile4 = tiled_partition(g, 2);
	// // if (tile4.thread_rank()==0) 
	// printf("Hello from tile4 rank %d: rank:%d\tthreadID:%d\tblockID:%d\n",tile4.thread_rank(),this_thread_block().thread_rank(),idx,level);

			
	__syncthreads();
	while((valid[idx*2 + (blockIdx.x)*N])==0 || (valid[idx*2 + (blockIdx.x)*N+1]) == 0);
		if(blockIdx.x == 0){
			in1 = network[idx*2];
			in2 = network[idx*2+1];
			if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x+1)*N] = in1;  
				network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
			}
			else{
				network[idx*2 + (blockIdx.x+1)*N] = in2;  
				network[idx*2 + (blockIdx.x+1)*N + 1] = in1;
			}
			valid[idx*2] = 0;  valid[idx*2 + 1] = 0;
			valid[idx*2 + (blockIdx.x+1)*N]=1; valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
			__syncthreads();
		}
		
		else {
			in1_index = LUT[idx*2 + (blockIdx.x-1)*N];
			in2_index = LUT[idx*2 + (blockIdx.x-1)*N + 1];
			in1 = network[in1_index+(blockIdx.x-1)*N];
			in2 = network[in2_index+(blockIdx.x-1)*N];
			if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x)*N] = in1;
				network[idx*2 + (blockIdx.x)*N + 1] = in2;
			}
			else{
				network[idx*2 + (blockIdx.x)*N] = in2;
				network[idx*2 + (blockIdx.x)*N + 1] = in1;  
			}
			valid[idx*2 + (blockIdx.x)*N] = 0; valid[idx*2 + 1 + (blockIdx.x)*N] = 0;
			valid[idx*2 + (blockIdx.x+1)*N]=1; valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;		
		} 
	

}



int main(int argc, char *argv[]){
	
	if (argc != 3){
		printf("Usage: %s <input.txt> <size>\n", argv[0]);
		return 1;
	}
	
	std::ifstream file(argv[1], std::ios::binary);
	if (!file) {
        printf("Could not open input file\n");
        return 1;
    }

	
	int N = atoi(argv[2]);
	int blockSize = N/2; 
	int numBlocks = 2*log2((double)N)-1; 
	int LUTsize = N*(log2((double)N)*2 - 2);
	
	char* network;
	cudaMallocManaged(&network,N*(numBlocks+1)*sizeof(char));
	memset(network,0,N*(numBlocks+1)*sizeof(char));
	file.read(network, N*sizeof(char));
	file.close();
	
	int* LUT;
	cudaMallocManaged(&LUT,LUTsize*sizeof(int));
	makeLUT(N,LUT);
	int mask = createMask(log2((double)N));
  
    int *valid;
	cudaMallocManaged(&valid,N*(numBlocks)*sizeof(int));
	memset(valid,0,N*(numBlocks+1)*sizeof(int)); 
	for(int i = 0; i < N; i++)
		valid[i] = 1;
	benes<<<numBlocks,blockSize>>>(N, network, LUT, valid, mask);
	cudaDeviceSynchronize();
	
	
	printf("The input is:");
	for (int i = 0; i < N; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", network[i]);
	}
	printf("\n");
	printf("The intermidiate layers are:\n");
	for (int i = N; i < N*(numBlocks-1); i++){
		if (i%N == 0) printf("\n");
		printf("%d ", network[i]);
		}
	printf("\n");
  
	for (int i = N*(numBlocks-1) ; i < N*(numBlocks); i++){
		if((mask & network[i]) != i % N){
			printf("ERROR in routing\n");
			return 1; 
		}
	}
	printf("Routing was successful!\nThe output is:\n");
	for (int i = N*(numBlocks-1); i < N*(numBlocks); i++){
		printf("%d ", network[i]);
	}
	printf("\n");
   
	cudaFree(valid);
	cudaFree(LUT);
	cudaFree(network);
}
 
 