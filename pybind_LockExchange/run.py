import numpy as np
import my_module

N = 10
arr = np.arange(N, dtype=np.float32)  # 初期条件 (0,1,2,...)

gpu_mem = my_module.GPUMemory(N)
# gpu_mem.set(arr)             # NumPy → GPU にコピー
my_module.compute_on_gpu(gpu_mem)  # CUDA カーネル呼び出し
result = gpu_mem.get()       # GPU → NumPy にコピー

print("Input :", arr)
print("Result:", result)
