#define __cplusplus
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "Fir_filter_cuda.cuh"
#include<stdio.h>
#include"device_functions.h"


__global__ void DoFIRFilt(const float* NUM, const int NUMLEN, const float* data, float* out, int datalen)
{
	extern __shared__ float tempsums[];
	const int tid = threadIdx.x;
	int did = blockIdx.x;
	for (size_t didx = did; didx < datalen; didx += gridDim.x)
	{
		float sum = 0;
		int i = tid;
		while (i < NUMLEN)
		{
			sum += NUM[i] * data[didx + i];
			i += blockDim.x;
		}
		tempsums[tid] = sum;
		__syncthreads();
		int numthreads = blockDim.x;
		//reduce
		if (numthreads >= 1024)
		{
			if (tid < 512)
				tempsums[tid] += tempsums[tid + 512];
			__syncthreads();
		}

		if (numthreads >= 512)
		{
			if (tid < 256)
				tempsums[tid] += tempsums[tid + 256];
			__syncthreads();
		}

		if (numthreads >= 256)
		{
			if (tid < 128)
				tempsums[tid] += tempsums[tid + 128];
			__syncthreads();
		}

		if (numthreads >= 128)
		{
			if (tid < 64)
				tempsums[tid] += tempsums[tid + 64];
			__syncthreads();
		}

		if (tid < 32)
		{
			volatile float* wssum = tempsums;
			if (numthreads >= 64)
				wssum[tid] += wssum[tid + 32];
			if (numthreads >= 32)
				wssum[tid] += wssum[tid + 16];
			if (numthreads >= 16)
				wssum[tid] += wssum[tid + 8];
			if (numthreads >= 8)
				wssum[tid] += wssum[tid + 4];
			if (numthreads >= 4)
				wssum[tid] += wssum[tid + 2];
			if (numthreads >= 2)
				wssum[tid] += wssum[tid + 1];

			if (tid == 0)
				out[didx] = wssum[0];
		}

		/*int j = blockDim.x/2;
		while (j > 0)
		{
			if(tid<j)
				tempsums[tid] += tempsums[tid + j];
			__syncthreads();
			j /= 2;
		}

		if (tid == 0)
			out[didx] = tempsums[0];*/
	}
}

extern "C" void fir_cuda(const float* NUM, int NUMLEN, float* data, int datalen, float* outputdata)
{
	int delay = (NUMLEN - 1) / 2;
	float* inputdata = new float[delay + datalen]();
	int i = 0;
	int inputlen = delay + datalen;
	while (i < datalen)
	{
		inputdata[i] = data[i];
		++i;
	}
	float* tempnum = new float[NUMLEN];
	i = 0;
	while (i < NUMLEN)
	{
		tempnum[i] = NUM[NUMLEN - 1 - i];
		++i;
	}
	cudaError_t cudaStatus;
	float* d_in = nullptr;
	float* d_out = nullptr;
	float* d_NUM = nullptr;
	cudaStatus = cudaMalloc((void**)& d_in, sizeof(float) * inputlen);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");

	}
	cudaStatus = cudaMalloc((void**)& d_out, sizeof(float) * datalen);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");

	}
	cudaStatus = cudaMalloc((void**)& d_NUM, sizeof(float) * NUMLEN);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");

	}
	cudaStatus = cudaMemcpy(d_in, inputdata, sizeof(float) * inputlen, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "copy failed!  Do you have a CUDA-capable GPU installed?");

	}
	cudaStatus = cudaMemcpy(d_NUM, tempnum, sizeof(float) * NUMLEN, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "copy failed!  Do you have a CUDA-capable GPU installed?");

	}

	int threadnum = NUMLEN;
	while ((threadnum & (threadnum - 1)) != 0)
		threadnum &= (threadnum - 1);
	if (threadnum > 2048)
		threadnum = 2048;
	int blocknum = (datalen < 100000 ? datalen : 100000);
	dim3 dimgrid(blocknum);
	dim3 blockdim(threadnum);
	size_t tnum = threadnum * sizeof(float);
	DoFIRFilt << <dimgrid, blockdim, tnum >> > (d_NUM, NUMLEN, d_in, d_out, datalen);


	cudaStatus = cudaMemcpy(outputdata, d_out, sizeof(float) * datalen, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");

	}
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaFree(d_in);
	cudaFree(d_out);
	delete[] inputdata;
	delete[] tempnum;
}