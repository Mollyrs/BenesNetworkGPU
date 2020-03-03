// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda.h>
#include <cuda_runtime.h>


//constucting 8x8 benes network
//four rows, 5 columns, 20 routers total

//put C:/Users/molly/Desktop/289Q/project/lut.cu
//nvcc -std=c++11 lut.cu
__global__
void benes(int N, int rows, int columns,  int* network){

}



int main(int argc, char *argv[]){
	int N = 32; //8x8 benes network

	int LUTsize = N*(log2((double)N)*2 - 2);

	int* LUT;
	cudaMallocManaged(&LUT,LUTsize*sizeof(int));
	
	int M = N;
	int even = 0;
	int odd = 1;
	
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
	
	
	//benes<<<numBlocks,blockSize>>>(N, 4, 5, network);
	//cudaDeviceSynchronize();
	
	for (int i = 0; i < LUTsize/2; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", LUT[i]);
	}
	printf("\n");
}