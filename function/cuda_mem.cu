#include "../all.hpp"

template<typename Typ>
void cuMallocCopy(Typ **d_val, const vector<Typ> &vec){
    cudaMalloc((void**)d_val, vec.size() * sizeof(Typ) ) ;
    cudaMemcpy(*d_val, vec.data(), vec.size() * sizeof(Typ), cudaMemcpyHostToDevice ) ;
}

template void cuMallocCopy(int**,    const vector<int>&);
template void cuMallocCopy(float**,  const vector<float>&);
template void cuMallocCopy(double**, const vector<double>&);