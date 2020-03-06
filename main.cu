//put C:/Users/molly/Desktop/289Q/project/main.cu
//nvcc -std=c++11 main.cu

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>

// includes, project
#include <cuda.h>
#include <cuda_runtime.h>


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


__global__
void benes(int N,  char* network, int* LUT){

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
	
	benes<<<numBlocks,blockSize>>>(N, network, LUT);
	cudaDeviceSynchronize();
	
	for (int i = 0; i < LUTsize; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", LUT[i]);
	}
	printf("\n");
	
	
	for (int i = 0; i < N*(numBlocks+1); i++){
		if (i%N == 0) printf("\n");
		printf("%d ", network[i]);
	}
	printf("\n");
}