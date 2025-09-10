#ifndef ALL_HPP      // まだ定義されていない場合に実行
#define ALL_HPP      // 定義する
using namespace std;

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <functional>
#include <future>
#include <type_traits>
#include <iostream>
#include <istream>
#include <memory>
#include <mutex>
#include <queue>
#include <stdio.h>
#include <sstream>
#include <string.h>
#include <string>
#include <thread>
#include <vector>

#include "itemIndex.hpp"

class Items{
    public:
        int nx, ny, nz, num_calc=0, num_IBMpoints=0, num_velocity ;
        int save_interval, total_count ;
        double dx, dt, ratiox, ratioy, PFthick, c, nu, sigma, tau, taus, gz, cx[27]={0}, cy[27]={0}, cz[27]={0}, weight[27] ;
        void setc(int alpha);
} ;

// basic.cu
template <typename Typ>
void output(const vector<Typ>& items, const vector<Typ>& posx, const vector<Typ>& posy, const vector<Typ>& posz, const vector<Typ>& delX, const vector<Typ>& delY,  const vector<Typ>& pressure, const vector<Typ>& velx, const vector<Typ>& vely, const vector<Typ>& velz, const vector<Typ>& sal, const vector<Typ>& phi, const vector<Typ>& rho, const vector<Typ>& Fx, const vector<Typ>& Fy, const vector<Typ>& Fz, int counter, int save_interval) ;
void input_items(Items &items, const char *fname) ;
template<typename Typ>
void reset_items_dt(vector<Typ> &items, Typ alpha) ;
__host__ __device__ void set_k_inv(int *k_inv, int num_velocity);
__host__ __device__ void set_cpm(int *cp, int *cm, int num_velocity);
__host__ __device__ void set_cpm(int *cp, int *cm, int *cip, int *cim, int num_velocity) ;
__host__ __device__ void set_cxypm(int *cxp, int *cxm, int *cyp, int *cym) ;

void set_constant_c(int num_velocity) ;
template <typename Typ>
void set_M(int num_velocity, vector<Typ> &M, vector<Typ> &S, vector<Typ> &M_inv, vector<Typ> &MM) ;
template<typename Typ>
void set_neibghor_wall(Items &items, int *lnum, vector<int> &neib, vector<Typ> &f, vector<Typ> &g, vector<Typ> &Fk, vector<Typ> &pressure, vector<Typ> &rho, vector<Typ> &phi, vector<Typ> &posx, vector<Typ> &posy, vector<Typ> &posz, vector<Typ> &velx, vector<Typ> &vely, vector<Typ> &velz) ;
template<typename Typ>
void set_neibghor_wall(Items &items, vector<int> &lnum, vector<int> &divx, vector<int> &divy, vector<int> &neib, vector<Typ> &f, vector<Typ> &g, vector<Typ> &Fk, vector<Typ> &pressure, vector<Typ> &rho, vector<Typ> &phi, vector<Typ> &posx, vector<Typ> &posy, vector<Typ> &posz, vector<Typ> &delX, vector<Typ> &delY, vector<Typ> &velx, vector<Typ> &vely, vector<Typ> &velz) ;
template<typename Typ>
void hydrostatic_pressure(Items &items, bool Boussi_flag, vector<int> &neib, vector<Typ> &pressure, vector<Typ> &rho, vector<Typ> &f, vector<Typ> &posz, Typ center_z) ;
template <typename Typ>
__global__ void set_wall_rho(Typ *items, int *neib, Typ *rho) ;
template<typename Typ>
__global__ void reset_pressure(Typ *items, Typ *pressure, Typ *f, Typ *velx, Typ *vely, Typ *velz, Typ alpha) ;
template<typename Typ>
__global__ void reset_wall_fin(Typ *items, Typ *fin, Typ *finy) ;
template<typename Typ>
__global__ void LES(Typ *items, int *neib, Typ *tau, Typ *taus, Typ *phi, Typ *rho, Typ muL, Typ muH, Typ *velx, Typ *vely, Typ *velz, Typ *posx, Typ *posy, Typ *posz) ;

// cuda_mem.cu
template<typename Typ>
void cuMallocCopy(Typ **d_val, const vector<Typ> &vec) ;

// collision.cu
template <typename Typ>
__global__ void equ_f(Typ *items, Typ *feq, Typ *pressure, Typ *u, Typ *v, Typ *w) ;
template <typename Typ>
__global__ void equ_g(Typ *items, Typ *feq, Typ *sal, Typ *u, Typ *v, Typ *w) ;
template <typename Typ>
__global__ void col_g_reg(Typ *items, Typ *taus, Typ *g, Typ *gtmp, Typ *geq, Typ *sal, Typ *velx, Typ *vely, Typ *velz) ;
template <typename Typ>
__global__ void col_f_SRT(Typ *items, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk) ;
template<typename Typ>
__global__ void col_f_SRT(Typ *items, Typ *tau, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk) ;
template <typename Typ>
__global__ void col_f_MRT(Typ *items, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk, Typ *M, Typ *M_inv, Typ *S) ;
template<typename Typ>
__global__ void col_f_MRT(Typ *items, Typ *tau, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk, Typ *M, Typ *M_inv, Typ *S, Typ *Mconst) ;
template<typename Typ>
__global__ void wall_function(Typ *items, Typ* delY, Typ* delX, int x, int y, int z, int wall_number, int *wall, Typ *velx, Typ *vely, Typ *velz, Typ *Fx, Typ *Fy, Typ *rho) ;
template <typename Typ>
__global__ void Force(Typ *items, bool Boussi_flag, int *neib, Typ *f, Typ *feq, Typ *tau, Typ *Fk, Typ *Fx, Typ *Fy, Typ *Fz, Typ *pressure, Typ *rho, Typ *sal, Typ *phi, Typ *velx, Typ *vely, Typ *velz, Typ *delX, Typ *delY, Typ *posx, Typ *posy, Typ *posz) ;
template<typename Typ>
__global__ void col_PF(Typ *items, int *neib, Typ *taus, Typ *g, Typ *gtmp, Typ *geq, Typ *phi, Typ *velx, Typ *vely, Typ *velz, Typ *phiold, Typ *uold, Typ *vold, Typ *wold, Typ *posx, Typ *posy, Typ *posz) ;
template<typename Typ>
__global__ void update_scalar(Typ *items, Typ *g, Typ *scalar) ;
template <typename Typ>
__global__ void update_rho(Typ *items, Typ rhoL, Typ rhoH, Typ *f, Typ *Fx, Typ *Fy, Typ *Fz, Typ *pressure, Typ *sal, Typ *phi, Typ *rho, Typ *u, Typ *v, Typ *w) ;
template<typename Typ>
__global__ void resetF(Typ *items, Typ *Fx, Typ *Fy, Typ *Fz, int num_F) ;

// ibm.cu
template<typename Typ>
void IB_csv(int loop, vector<Typ>& items, vector<Typ>& pos, vector<Typ>& velw, vector<Typ>& Gw) ;
__host__ __device__ void set_quaternionS(int IB_index, float q0, float q1, float q2, float q3, float *quaS) ;
template<typename Typ>
void set_quaternionS(int IB_index, Typ q0, Typ q1, Typ q2, Typ q3, vector<Typ>& quaS) ;
template<typename Typ>
__global__ void SPM(Typ *items, Typ Radius, Typ *posB, Typ *f, Typ *ftmp, Typ *tau, Typ *posx, Typ *posy, Typ *Fx, Typ *Fy, Typ *Fz, Typ *velx, Typ *vely, Typ *velz, Typ *velw) ;
template<typename Typ>
__global__ void SPM_ellipse(Typ *items, Typ Rada, Typ Radb, Typ *quaS, Typ *posB, Typ *f, Typ *tau, Typ *posx, Typ *posy, Typ *posz, Typ *velx, Typ *vely, Typ *velz, Typ *velB, Typ *angleVB) ;
template<typename Typ>
__global__ void get_IBMGw2(Typ *items, int *lattice_id, int *neib, Typ *f, Typ *tau, Typ *posx, Typ *posy, Typ *posz, Typ *posw, Typ *posB, Typ *nBvec, Typ *velx, Typ *vely, Typ *velz, Typ *velw, Typ *Fx, Typ *Fy, Typ *Fz, Typ *Gw, Typ rhof) ;
__global__ void IB_directForcing(float *items, int *lattice_id, int *neib, float *posx, float *posy, float *posz, float *posw,float *velx, float *vely, float *velz, float *velw, float *Fx, float *Fy, float *Fz, float *Gw) ;
template<typename Typ>
__global__ void update_velIBM(Typ *items, int *lattice_id, Typ *f, Typ *ftmp, Typ *pressure, Typ *tau, Typ *velx, Typ *vely, Typ *velz, Typ *velx_old, Typ *vely_old, Typ *velz_old, Typ *Fx, Typ *Fy, Typ *Fz) ;
__global__ void update_IBbody(float *items, int IB_index, float *massB, float *densB, float *inertia, float *FB, float *posB, float *Torque, float *velB, float *quat, float *quaS, float *angleVB, float *posw, float *Gw, float *quatold, float rhof) ;
__global__ void update_IBpoint(float *items, int IB_index, float *posB, float *velB, float *angleVB, float *quaS, float *posw, float *oposw, float *nBvec, float *onBvec, float *velw) ;
__global__ void search_IBlattice(float *items, int IB_index, int *lattice_id, int *neib, float *posx, float *posy, float *posz, float *posw) ;



// adv.cu
template <typename Typ>
__global__ void set_wall_f(Typ *items, int *neib, Typ *f) ;
template <typename Typ>
__global__ void set_out(Typ *items, int *neib, Typ *ftmp, Typ *fout, Typ *fin, Typ *finy, Typ *delX, Typ *delY) ;
template<typename Typ>
__global__ void set_fin(Typ *items, int *neib, int *nextK, Typ *fin, Typ *finy) ;
template<typename Typ>
__global__ void reset_f(Typ *items, Typ *ftmp, Typ *fin, Typ *f) ;
template<typename Typ>
__global__ void set_wall_f_inv(Typ *items, int *neib, int *nextB, int *nextK, Typ *f) ;
template<typename Typ>
__device__ void  Lagrange_interpolation(Typ *coeff, Typ *posl_x, Typ x_alpha, int num_point) ;
template<typename Typ>
void IP_process(Typ *d_items, int numblocks, int blockSize, int *d_neib, Typ *d_f, Typ *d_feq, Typ *d_ftmp, Typ *d_fout, int *d_nextB, int *d_nextK, Typ *d_posx, Typ *d_posy, Typ *d_delX, Typ *d_delY, int slip) ;
template <typename Typ>
__global__ void set_tmp(Typ *items, int *neib, Typ *f, Typ *ftmp, Typ *fout, Typ *fin) ;

// set_bound
template<typename Typ>
void set_walln(const vector<Typ>& items, const vector<int>& neib, vector<int>& wall1, vector<int>& wall2, vector<int>& wall3, vector<int>& wall4, vector<int>& wall5, vector<int>& wall6) ;
template <typename Typ>
void set_bound2D(vector<Typ>& items, int num_block, const vector<int>& neib, vector<int>& nextK, vector<int>& nextB) ;
template<typename Typ>
void set_bound3D(vector<Typ>& items, int num_block, const vector<int>& neib, vector<int>& nextK, vector<int>& nextB) ;
template <typename Typ>
__global__ void set_wall_propagation(Typ *items, int boudrary_condition, int *neib, int *nextK, int *nextB, Typ *ftmp) ;

// propagation.cu
template <typename Typ>
__global__ void propagation(Typ *items, int *neib, Typ *f, Typ *ftmp) ;


#endif