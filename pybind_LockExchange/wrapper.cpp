#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cuda_runtime.h>

namespace py = pybind11;

extern "C" void launch_add_one(float* , int ) ;


// GPU メモリを保持するクラス
class GPUMemory {
    float* d_data;
    size_t N;

public:
    GPUMemory(size_t n) : N(n) {
        cudaMalloc(&d_data, N * sizeof(float));
    }

    ~GPUMemory() {
        cudaFree(d_data);
    }

    // ホストからGPUにコピー
    void set(py::array_t<float> arr) {
        py::buffer_info buf = arr.request();
        if ((size_t)buf.size != N) throw std::runtime_error("Size mismatch");
        cudaMemcpy(d_data, buf.ptr, N * sizeof(float), cudaMemcpyHostToDevice);
    }

    // GPUからホストにコピー
    py::array_t<float> get() {
        py::array_t<float> result(N);
        py::buffer_info buf = result.request();
        cudaMemcpy(buf.ptr, d_data, N * sizeof(float), cudaMemcpyDeviceToHost);
        return result;
    }

    float* data() { return d_data; }
    size_t size() { return N; }
};

// CUDAカーネル呼び出し（kernel.cu で定義）
extern void launch_add_one(float* d_data, int N);

void cumallcpy()

void compute_on_gpu(GPUMemory& mem) {
    launch_add_one(mem.data(), mem.size());
}

// Python モジュール定義
PYBIND11_MODULE(my_module, m) {
    py::class_<GPUMemory>(m, "GPUMemory")
        .def(py::init<size_t>())
        .def("set", &GPUMemory::set)
        .def("get", &GPUMemory::get)
        .def("size", &GPUMemory::size);

    m.def("compute_on_gpu", &compute_on_gpu);
}
