#include <cuda_runtime.h>

// 簡単なカーネル: 配列の各要素に +1
__global__ void add_one(float* d_data, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        d_data[idx] += 5.0f;
    }
}

// GPU 計算を呼び出す関数
extern "C" void launch_add_one(float* d_data, int N) {
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    add_one<<<numBlocks, blockSize>>>(d_data, N);
    cudaDeviceSynchronize();
}
