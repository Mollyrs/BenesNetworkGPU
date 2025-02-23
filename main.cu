//put C:/Users/molly/Desktop/289Q/project/main.cu
//nvcc -std=c++11 main.cu

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <cooperative_groups.h>
#include <cooperative_groups.h>
// includes, project
#include <cuda.h>
#include <cuda_runtime.h>

using namespace cooperative_groups;
namespace cg = cooperative_groups;

// #define FILESIZE_CHAR 1048576
#define FILESIZE_CHAR 1048576
#define FILESIZE_INT FILESIZE_CHAR/4


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
void benes(int N, int block, char* network, int* LUT, volatile int* valid, int mask, int* data, char* output){
	int idx = threadIdx.x;
	int in1, in2, in1_index, in2_index;
	int readOffset=0;
	int fileSize = FILESIZE_INT/2;
	int readOffsetSecondNet=fileSize;
	thread_group g = tiled_partition(this_thread_block(), 2); //stops working after 32?
		if(blockIdx.x == 0){
			while(readOffset < fileSize){
				in1 = data[idx*2 + readOffset];
				in2 = data[idx*2+1 + readOffset];
				readOffset+=N;
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x+1)*N] = in1;  
				network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
				
				}
				else{
					network[idx*2 + (blockIdx.x+1)*N] = in2;  
					network[idx*2 + (blockIdx.x+1)*N + 1] = in1;
				}
				g.sync();
				// __syncthreads();
				valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
				
			}
		}
		
		else if ( blockIdx.x < block) {
			while(readOffset < fileSize){
				while((valid[idx + (blockIdx.x)*(N/2)])==0);
				in1_index = LUT[idx*2 + (blockIdx.x-1)*N];
				in2_index = LUT[idx*2 + (blockIdx.x-1)*N + 1];
				in1 = network[in1_index+(blockIdx.x)*N];
				in2 = network[in2_index+(blockIdx.x)*N];
				valid[idx + (blockIdx.x)*(N/2)] = 0;// valid[idx*2 + 1 + (blockIdx.x)*N] = 0;
			

				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x+1)*N] = in1;
				network[idx*2 + (blockIdx.x+1)*N + 1] = in2;

				}
				else{
					network[idx*2 + (blockIdx.x+1)*N] = in2;
					network[idx*2 + (blockIdx.x+1)*N + 1] = in1;  
				}
				
				if (blockIdx.x != gridDim.x - 1 && blockIdx.x != block-1){
					valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
					g.sync();
					// __syncthreads();
				}
				else {
					output[idx*2 + readOffset] = network[idx*2 + (blockIdx.x+1)*N];
					output[idx*2+1 + readOffset] = network[idx*2 + (blockIdx.x+1)*N + 1];
				}
				readOffset += N;
			}
		} 



		else if(blockIdx.x == block){
			while(readOffsetSecondNet < FILESIZE_INT){
				in1 = data[idx*2 + readOffsetSecondNet];
				in2 = data[idx*2+1 + readOffsetSecondNet];
				readOffsetSecondNet+=N;
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				if ((in1 & mask) < (in2 & mask)){
					network[idx*2 + (blockIdx.x+1)*N] = in1;  
					network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
				}
				else{
					network[idx*2 + (blockIdx.x+1)*N] = in2;  
					network[idx*2 + (blockIdx.x+1)*N + 1] = in1;
				}

				
				valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
				// __syncthreads();
				g.sync();
			}
		}
		
		else{
			while(readOffsetSecondNet < FILESIZE_INT){
				// printf("waiting for previous block %d to produce\n", blockIdx.x - 1);
				while((valid[idx + (blockIdx.x)*(N/2)])==0);
				
				// printf("waiting for previous block %d to produce\n", blockIdx.x - 1);
				in1_index = LUT[idx*2 + ((blockIdx.x%block)-1)*N];
				in2_index = LUT[idx*2 + ((blockIdx.x%block)-1)*N + 1];
				in1 = network[in1_index+(blockIdx.x)*N];
				in2 = network[in2_index+(blockIdx.x)*N];
				
				// printf("Block %d thread %d consumed %d %d\n", blockIdx.x,threadIdx.x, in1, in2);
				valid[idx + (blockIdx.x)*(N/2)] = 0; //valid[idx*2 + 1 + (blockIdx.x)*N] = 0;
			
				//printf("waiting for next block %d to consume\n", blockIdx.x + 1);
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				if ((in1 & mask) < (in2 & mask)){
					network[idx*2 + (blockIdx.x+1)*N] = in1;
					network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
					// printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				}
				else{
					network[idx*2 + (blockIdx.x+1)*N] = in2;
					network[idx*2 + (blockIdx.x+1)*N + 1] = in1;  
				}
				//printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				if (blockIdx.x != gridDim.x - 1){
					valid[idx + (blockIdx.x+1)*(N/2)]=1; //valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
					// __syncthreads();
					g.sync();
					//printf("valid:%d index:%d\n",valid[idx + (blockIdx.x+1)*N],idx + (blockIdx.x+1)*N);
				}
				else {
					output[idx*2 + readOffsetSecondNet] = network[idx*2 + (blockIdx.x+1)*N];
					output[idx*2+1 + readOffsetSecondNet] = network[idx*2 + (blockIdx.x+1)*N + 1];
				}
				readOffsetSecondNet += N;
			}
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
	if (FILESIZE_INT<N)
		N = FILESIZE_INT;
	int blockSize = N/2; 
	int blocks = 2*log2((double)N)-1; 
	int b = 2*log2((double)N)-1;
	int LUTsize = N*(log2((double)N)*2 - 2);
	int numBlocks;

	if (FILESIZE_INT <= N)
		numBlocks = blocks;
	else
		numBlocks = 2*blocks;

	char* network;
	cudaMallocManaged(&network,N*(numBlocks+1)*sizeof(char));
	memset(network,0,N*(numBlocks+1)*sizeof(char));

	
	int* LUT;
	cudaMallocManaged(&LUT,LUTsize*sizeof(int));
	makeLUT(N,LUT);
	int mask = createMask(log2((double)N));
  
    int *valid;
	cudaMallocManaged(&valid,(N/2)*(numBlocks)*sizeof(int));
	memset(valid,0,(N/2)*(numBlocks+1)*sizeof(int)); 
	for(int i = 0; i < N/2; i++)
		valid[i] = 1;
	
	char* data;
	cudaMallocManaged(&data,FILESIZE_CHAR*sizeof(char));
	memset(data,0,FILESIZE_CHAR*sizeof(char));
	file.read(data, FILESIZE_CHAR*sizeof(char));
	file.close();
	
	int* idata;
	cudaMallocManaged(&idata,FILESIZE_CHAR*sizeof(char));
	memcpy(idata, data, FILESIZE_CHAR*sizeof(char));
	
	char* output;
	cudaMallocManaged(&output,FILESIZE_CHAR*sizeof(char));
	memset(output,0,FILESIZE_CHAR*sizeof(char));
	
	benes<<<numBlocks,blockSize>>>(N, blocks, network, LUT, valid, mask, idata, output);
	cudaDeviceSynchronize();

	// printf("The input is:");
	// for (int i = 0; i < FILESIZE_INT; i++){
		// if (i%N == 0) printf("\n");
		// printf("%d ", idata[i]);
	// }
	// printf("\n\n");

  
	for (int i = 0; i < FILESIZE_INT-1; i++){
		if ((i%N != N-1) && (output[i+1]!=0)) {
				if((mask & output[i+1]) < (mask & output[i])){
						printf("ERROR in routing at output %d %d %d\n",i ,mask & output[i+1],mask &output[i] );
						return 1;
				}
		}
}
printf("Routing was successful!\n");
   
	cudaFree(valid);
	cudaFree(LUT);
	cudaFree(network);
	cudaFree(data);
	cudaFree(idata);
	cudaFree(output);
}
 
 
