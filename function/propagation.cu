#include "../all.hpp"
template<typename Typ>
__global__ void propagation(Typ *items, int *neib, Typ *f, Typ *ftmp){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*int(items[IDX_Q]) ;
    if(id_rho<items[IDX_num_calc]){
        int k, k_inv[27]={0} ; set_k_inv(k_inv,(int)items[IDX_Q]) ;
        for(k=0;k<items[IDX_Q];k++){
            f[id_f+k] = ftmp[neib[id_f+k_inv[k]]*(int)items[IDX_Q]+k] ;
        }
    }// */
}

template __global__ void propagation<float>(float*, int*, float*, float*) ;

template __global__ void propagation<double>(double*, int*, double*, double*) ;