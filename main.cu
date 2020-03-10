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

#define FILESIZE 256



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
void benes(int N, int block, char* network, int* LUT, volatile int* valid, int mask, char* data, char* output){
	int idx = threadIdx.x;
	int in1, in2, in1_index, in2_index;
	int readOffset=0;
	int fileSize = FILESIZE/2;
	int readOffsetSecondNet=fileSize;
	thread_group g = tiled_partition(this_thread_block(), 2);
		if(blockIdx.x == 0){
			while(readOffset < fileSize){
				in1 = data[idx*2 + readOffset];
				in2 = data[idx*2+1 + readOffset];
				readOffset+=N;
				//printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				//printf("waiting for next block %d to consume\n", blockIdx.x + 1);
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				// if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x+1)*N] = in1;  
				network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
				
				// }
				// else{
				// 	network[idx*2 + (blockIdx.x+1)*N] = in2;  
				// 	network[idx*2 + (blockIdx.x+1)*N + 1] = in1;
				// }
				g.sync();
				valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
				
			}
		}
		
		else if ( blockIdx.x < block) {
			while(readOffset < fileSize){
				//printf("waiting for previous block %d to produce\n", blockIdx.x - 1);
				while((valid[idx + (blockIdx.x)*(N/2)])==0);
				in1_index = LUT[idx*2 + (blockIdx.x-1)*N];
				in2_index = LUT[idx*2 + (blockIdx.x-1)*N + 1];
				in1 = network[in1_index+(blockIdx.x)*N];
				in2 = network[in2_index+(blockIdx.x)*N];
				//printf("Block %d consumed %d %d\n", blockIdx.x, in1, in2);
				valid[idx + (blockIdx.x)*(N/2)] = 0;// valid[idx*2 + 1 + (blockIdx.x)*N] = 0;
			
				//printf("waiting for next block %d to consume\n", blockIdx.x + 1);
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				// if ((in1 & mask) < (in2 & mask)){
				network[idx*2 + (blockIdx.x+1)*N] = in1;
				network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
				g.sync();
				// }
				// else{
				// 	network[idx*2 + (blockIdx.x+1)*N] = in2;
				// 	network[idx*2 + (blockIdx.x+1)*N + 1] = in1;  
				// }
				
				if (blockIdx.x != gridDim.x - 1 && blockIdx.x != block-1){
					g.sync();
					valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
				}
				else {
					output[idx*2 + readOffset] = network[idx*2 + (blockIdx.x+1)*N];
					output[idx*2+1 + readOffset] = network[idx*2 + (blockIdx.x+1)*N + 1];
				}
				// printf("Block %d produced %d %d\n", gridDim.x,output[idx*2 + readOffset], output[idx*2+1 + readOffset]);
				readOffset += N;
			}
		} 



		else if(blockIdx.x == block){
			while(readOffsetSecondNet < FILESIZE){
				in1 = data[idx*2 + readOffsetSecondNet];
				in2 = data[idx*2+1 + readOffsetSecondNet];
				readOffsetSecondNet+=N;
				// printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				//printf("waiting for next block %d to consume\n", blockIdx.x + 1);
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				// if ((in1 & mask) < (in2 & mask)){
					network[idx*2 + (blockIdx.x+1)*N] = in1;  
					network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
				// }
				// else{
				// 	network[idx*2 + (blockIdx.x+1)*N] = in2;  
				// 	network[idx*2 + (blockIdx.x+1)*N + 1] = in1;
				// }
				__syncthreads();
				// printf("Block %d produced %d %d\n", blockIdx.x, network[idx*2 + (blockIdx.x+1)*N],network[idx*2 + (blockIdx.x+1)*N+1]);
				valid[idx + (blockIdx.x+1)*(N/2)]=1;// valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
				__syncthreads();
			}
		}
		
		else{
			while(readOffsetSecondNet < FILESIZE){
				// printf("waiting for previous block %d to produce\n", blockIdx.x - 1);
				while((valid[idx + (blockIdx.x)*(N/2)])==0);
				__syncthreads();
				
				// printf("waiting for previous block %d to produce\n", blockIdx.x - 1);
				in1_index = LUT[idx*2 + ((blockIdx.x%block)-1)*N];
				in2_index = LUT[idx*2 + ((blockIdx.x%block)-1)*N + 1];
				in1 = network[in1_index+(blockIdx.x)*N];
				in2 = network[in2_index+(blockIdx.x)*N];
				
				// printf("Block %d thread %d consumed %d %d\n", blockIdx.x,threadIdx.x, in1, in2);
				valid[idx + (blockIdx.x)*(N/2)] = 0; //valid[idx*2 + 1 + (blockIdx.x)*N] = 0;
			
				//printf("waiting for next block %d to consume\n", blockIdx.x + 1);
				while((valid[idx + (blockIdx.x+1)*(N/2)])==1);
				// if ((in1 & mask) < (in2 & mask)){
					network[idx*2 + (blockIdx.x+1)*N] = in1;
					network[idx*2 + (blockIdx.x+1)*N + 1] = in2;
					// printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				// }
				// else{
				// 	network[idx*2 + (blockIdx.x+1)*N] = in2;
				// 	network[idx*2 + (blockIdx.x+1)*N + 1] = in1;  
				// }
				//printf("Block %d produced %d %d\n", blockIdx.x, in1, in2);
				if (blockIdx.x != gridDim.x - 1){
					valid[idx + (blockIdx.x+1)*(N/2)]=1; //valid[idx*2 + 1 + (blockIdx.x+1)*N]=1;
					printf("valid:%d index:%d\n",valid[idx + (blockIdx.x+1)*N],idx + (blockIdx.x+1)*N);
				}
				else {
					output[idx*2 + readOffsetSecondNet] = network[idx*2 + (blockIdx.x+1)*N];
					output[idx*2+1 + readOffsetSecondNet] = network[idx*2 + (blockIdx.x+1)*N + 1];
				}
				// printf("HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");
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
	if (FILESIZE<N)
		N = FILESIZE;
	int blockSize = N/2; 
	int blocks = 2*log2((double)N)-1; 
	int b = 2*log2((double)N)-1;
	int LUTsize = N*(log2((double)N)*2 - 2);
	int numBlocks;

	if (FILESIZE <= N)
		numBlocks = blocks;
	else
		numBlocks = 2*blocks;

	char* network;
	cudaMallocManaged(&network,N*(numBlocks+1)*sizeof(char));
	memset(network,0,N*(numBlocks+1)*sizeof(char));
	//file.read(network, N*sizeof(char));
	//file.close();
	
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
	cudaMallocManaged(&data,FILESIZE*sizeof(char));
	memset(data,0,FILESIZE*sizeof(char));
	file.read(data, FILESIZE*sizeof(char));
	file.close();
	
	char* output;
	cudaMallocManaged(&output,FILESIZE*sizeof(char));
	memset(output,0,FILESIZE*sizeof(char));

	
	
	benes<<<numBlocks,blockSize>>>(N, blocks, network, LUT, valid, mask, data, output);
	cudaDeviceSynchronize();
	
	
	
	
	printf("The input is:");
	for (int i = 0; i < FILESIZE; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", data[i]);
	}
	printf("\n\n");

  
	printf("The output is:");
	for (int i = 0; i < FILESIZE; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", output[i]);
	}
	printf("\n");
   
	cudaFree(valid);
	cudaFree(LUT);
	cudaFree(network);
	cudaFree(data);
	cudaFree(output);
}
 
 
