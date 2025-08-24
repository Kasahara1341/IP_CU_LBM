#include "../all.hpp"
// num_block : 計算格子の総数
// neib  : 隣接格子の番号
// nextK : 取ってくるkの番号
// nextB : 計算格子から見た取ってくる格子の方向
// nextKとnextBのiは壁格子としてのインデックス(0から始まり壁格子の総数で終わる)

template<typename Typ>
void set_walln(const vector<Typ>& items, const vector<int>& neib, vector<int>& wall1, vector<int>& wall2, vector<int>& wall3, vector<int>& wall4, vector<int>& wall5, vector<int>& wall6){
    int i ;
    int dir1=1, dir2, dir3=3, dir4, dir5, dir6 ;
    if      ((int)items[IDX_Q]==9 ){ dir2=0 ; dir4=0 ; dir5=2 ; dir6=4 ;printf("num_velocity = %d\n",(int)items[IDX_Q]);}
    else if ((int)items[IDX_Q]==27){ dir2=2 ; dir4=4 ; dir5=5 ; dir6=6 ;printf("num_velocity = %d\n",(int)items[IDX_Q]);}
    for(i=0;i<items[IDX_num_calc];i++){
        if(neib[i*(int)items[IDX_Q]+dir1]>=(int)items[IDX_num_calc]){ wall1.push_back(i) ;}
        if(neib[i*(int)items[IDX_Q]+dir2]>=(int)items[IDX_num_calc]){ wall2.push_back(i) ;}
        if(neib[i*(int)items[IDX_Q]+dir3]>=(int)items[IDX_num_calc]){ wall3.push_back(i) ;}
        if(neib[i*(int)items[IDX_Q]+dir4]>=(int)items[IDX_num_calc]){ wall4.push_back(i) ;}
        if(neib[i*(int)items[IDX_Q]+dir5]>=(int)items[IDX_num_calc]){ wall5.push_back(i) ;}
        if(neib[i*(int)items[IDX_Q]+dir6]>=(int)items[IDX_num_calc]){ wall6.push_back(i) ;}

    }
}

template<typename Typ>
void set_bound2D(vector<Typ>& items, int num_block, const vector<int>& neib, vector<int>& nextK, vector<int>& nextB){
    int i, k, cxx, czz, cx[9], cz[9] ; 
    int nei[2][2][9]={0}, kinv[2][2][9], k_inv[27]={0} ; set_k_inv(k_inv,(int)items[IDX_Q]) ;
    for(k=0;k<9;k++){
        cx[k]=items[IDX_cx(k)]/items[IDX_c] ; cz[k]=items[IDX_cz(k)]/items[IDX_c] ;
        nei[abs(cx[k])][abs(cz[k])][k]=0 ; kinv[abs(cx[k])][abs(cz[k])][k]=k_inv[k] ; kinv[0][0][k]=k_inv[k] ;
    }

    nei[1][0][5] = 3 ; nei[0][1][5] = 4 ; nei[1][0][6] = 1 ; nei[0][1][6] = 4 ;
    nei[1][0][7] = 1 ; nei[0][1][7] = 2 ; nei[1][0][8] = 3 ; nei[0][1][8] = 2 ;
    
    kinv[1][0][5] = 8 ; kinv[0][1][5] = 6 ; kinv[1][0][6] = 7 ; kinv[0][1][6] = 5 ;
    kinv[1][0][7] = 6 ; kinv[0][1][7] = 8 ; kinv[1][0][8] = 5 ; kinv[0][1][8] = 7 ;

    for(i=0;i<items[IDX_num_wall];i++){
        int calc_id=neib[(num_block+i)*(int)items[IDX_Q]]*(int)items[IDX_Q] ;
        nextB.push_back(0) ;
        nextK.push_back(0) ;
        for(k=1;k<9;k++){
            cxx=cx[k] ; czz=cz[k] ;
            if(neib[calc_id+k_inv[k]]>=num_block){
                if((neib[calc_id+1]>=num_block && cxx==-1) || (neib[calc_id+3]>=num_block && cxx==1)){
                    cxx=0 ;
                }
                if((neib[calc_id+2]>=num_block && czz==-1) || (neib[calc_id+4]>=num_block && czz==1)){
                    czz=0 ;
                }
            }
            nextB.push_back(nei[abs(cxx)][abs(czz)][k]) ;
            nextK.push_back(kinv[abs(cxx)][abs(czz)][k]) ;
        }
    }
}

template<typename Typ>
void set_bound3D(vector<Typ>& items, int num_block, const vector<int>& neib, vector<int>& nextK, vector<int>& nextB){
    int i, k, cxx, cyy, czz, cx[27], cy[27], cz[27] ; 
    int nei[2][2][2][27]={0}, kinv[2][2][2][27], k_inv[27] ; set_k_inv(k_inv,(int)items[IDX_Q]) ;

    for(k=0;k<27;k++){
        cx[k]=items[IDX_cx(k)]/items[IDX_c] ; cy[k]=items[IDX_cy(k)]/items[IDX_c] ; cz[k]=items[IDX_cz(k)]/items[IDX_c] ;
        nei[abs(cx[k])][abs(cy[k])][abs(cz[k])][k]=0 ; kinv[abs(cx[k])][abs(cy[k])][abs(cz[k])][k]=k_inv[k] ; kinv[0][0][0][k]=k_inv[k] ;
    }

    nei[1][0][0][7 ] = 3 ; nei[0][1][0][7 ] = 4 ; nei[1][0][0][8 ] = 3 ; nei[0][1][0][8 ] = 2 ;
    nei[1][0][0][9 ] = 1 ; nei[0][1][0][9 ] = 4 ; nei[1][0][0][10] = 1 ; nei[0][1][0][10] = 2 ;
    nei[1][0][0][11] = 3 ; nei[0][0][1][11] = 6 ; nei[1][0][0][12] = 3 ; nei[0][0][1][12] = 5 ;
    nei[1][0][0][13] = 1 ; nei[0][0][1][13] = 6 ; nei[1][0][0][14] = 1 ; nei[0][0][1][14] = 5 ;
    nei[0][1][0][15] = 4 ; nei[0][0][1][15] = 6 ; nei[0][1][0][16] = 4 ; nei[0][0][1][16] = 5 ;
    nei[0][1][0][17] = 2 ; nei[0][0][1][17] = 6 ; nei[0][1][0][18] = 2 ; nei[0][0][1][18] = 5 ;
    nei[1][0][0][19] = 3 ; nei[0][1][0][19] = 4 ; nei[0][0][1][19] = 6 ; nei[1][1][0][19] = 10 ; nei[1][0][1][19] = 14 ; nei[0][1][1][19] = 18 ; 
    nei[1][0][0][20] = 3 ; nei[0][1][0][20] = 4 ; nei[0][0][1][20] = 5 ; nei[1][1][0][20] = 10 ; nei[1][0][1][20] = 13 ; nei[0][1][1][20] = 17 ;
    nei[1][0][0][21] = 3 ; nei[0][1][0][21] = 2 ; nei[0][0][1][21] = 6 ; nei[1][1][0][21] = 9  ; nei[1][0][1][21] = 14 ; nei[0][1][1][21] = 16 ;
    nei[1][0][0][22] = 3 ; nei[0][1][0][22] = 2 ; nei[0][0][1][22] = 5 ; nei[1][1][0][22] = 9  ; nei[1][0][1][22] = 13 ; nei[0][1][1][22] = 15 ;
    nei[1][0][0][23] = 1 ; nei[0][1][0][23] = 4 ; nei[0][0][1][23] = 6 ; nei[1][1][0][23] = 8  ; nei[1][0][1][23] = 12 ; nei[0][1][1][23] = 18 ;
    nei[1][0][0][24] = 1 ; nei[0][1][0][24] = 4 ; nei[0][0][1][24] = 5 ; nei[1][1][0][24] = 8  ; nei[1][0][1][24] = 11 ; nei[0][1][1][24] = 17 ;
    nei[1][0][0][25] = 1 ; nei[0][1][0][25] = 2 ; nei[0][0][1][25] = 6 ; nei[1][1][0][25] = 7  ; nei[1][0][1][25] = 12 ; nei[0][1][1][25] = 16 ;
    nei[1][0][0][26] = 1 ; nei[0][1][0][26] = 2 ; nei[0][0][1][26] = 5 ; nei[1][1][0][26] = 7  ; nei[1][0][1][26] = 11 ; nei[0][1][1][26] = 15 ;

    kinv[1][0][0][7 ] = 8  ; kinv[0][1][0][7 ] = 9  ; kinv[1][0][0][8 ] = 7  ; kinv[0][1][0][8 ] = 10 ;
    kinv[1][0][0][9 ] = 10 ; kinv[0][1][0][9 ] = 7  ; kinv[1][0][0][10] = 9  ; kinv[0][1][0][10] = 8  ;
    kinv[1][0][0][11] = 12 ; kinv[0][0][1][11] = 13 ; kinv[1][0][0][12] = 11 ; kinv[0][0][1][12] = 14 ;
    kinv[1][0][0][13] = 14 ; kinv[0][0][1][13] = 11 ; kinv[1][0][0][14] = 13 ; kinv[0][0][1][14] = 12 ;
    kinv[0][1][0][15] = 16 ; kinv[0][0][1][15] = 17 ; kinv[0][1][0][16] = 15 ; kinv[0][0][1][16] = 18 ;
    kinv[0][1][0][17] = 18 ; kinv[0][0][1][17] = 15 ; kinv[0][1][0][18] = 17 ; kinv[0][0][1][18] = 16 ;
    kinv[1][0][0][19] = 22 ; kinv[0][1][0][19] = 24 ; kinv[0][0][1][19] = 25 ; kinv[1][1][0][19] = 20 ; kinv[1][0][1][19] = 21 ; kinv[0][1][1][19] = 23 ; 
    kinv[1][0][0][20] = 21 ; kinv[0][1][0][20] = 23 ; kinv[0][0][1][20] = 26 ; kinv[1][1][0][20] = 19 ; kinv[1][0][1][20] = 22 ; kinv[0][1][1][20] = 24 ;
    kinv[1][0][0][21] = 20 ; kinv[0][1][0][21] = 26 ; kinv[0][0][1][21] = 23 ; kinv[1][1][0][21] = 22 ; kinv[1][0][1][21] = 19 ; kinv[0][1][1][21] = 25 ;
    kinv[1][0][0][22] = 19 ; kinv[0][1][0][22] = 25 ; kinv[0][0][1][22] = 24 ; kinv[1][1][0][22] = 21 ; kinv[1][0][1][22] = 20 ; kinv[0][1][1][22] = 26 ;
    kinv[1][0][0][23] = 26 ; kinv[0][1][0][23] = 20 ; kinv[0][0][1][23] = 21 ; kinv[1][1][0][23] = 24 ; kinv[1][0][1][23] = 25 ; kinv[0][1][1][23] = 19 ;
    kinv[1][0][0][24] = 25 ; kinv[0][1][0][24] = 19 ; kinv[0][0][1][24] = 22 ; kinv[1][1][0][24] = 23 ; kinv[1][0][1][24] = 26 ; kinv[0][1][1][24] = 20 ;
    kinv[1][0][0][25] = 24 ; kinv[0][1][0][25] = 22 ; kinv[0][0][1][25] = 19 ; kinv[1][1][0][25] = 26 ; kinv[1][0][1][25] = 23 ; kinv[0][1][1][25] = 21 ;
    kinv[1][0][0][26] = 23 ; kinv[0][1][0][26] = 21 ; kinv[0][0][1][26] = 20 ; kinv[1][1][0][26] = 25 ; kinv[1][0][1][26] = 24 ; kinv[0][1][1][26] = 22 ;    

    for(i=0;i<items[IDX_num_wall];i++){
        int calc_id=neib[(num_block+i)*(int)items[IDX_Q]]*(int)items[IDX_Q] ;
        nextB.push_back(0) ;
        nextK.push_back(0) ;
        for(k=1;k<items[IDX_Q];k++){
            cxx=cx[k] ; cyy=cy[k] ; czz=cz[k] ;
            if(neib[calc_id+k_inv[k]]>=num_block){
                if((neib[calc_id+1]>=num_block && cxx==-1) || (neib[calc_id+3]>=num_block && cxx==1)){
                    cxx=0 ;
                }
                if((neib[calc_id+2]>=num_block && cyy==-1) || (neib[calc_id+4]>=num_block && cyy==1)){
                    cyy=0 ;
                }
                if((neib[calc_id+5]>=num_block && czz==-1) || (neib[calc_id+6]>=num_block && czz==1)){
                    czz=0 ;
                }

                if((neib[calc_id+11]<num_block && cx[k]==-1 && cz[k]==-1) || (neib[calc_id+14]<num_block && cx[k]==1 && cz[k]==1)){
                    cxx=cx[k] ; czz=cz[k] ;
                }
                if((neib[calc_id+12]<num_block && cx[k]==-1 && cz[k]==1) || (neib[calc_id+13]<num_block && cx[k]==1 && cz[k]==-1)){
                    cxx=cx[k] ; czz=cz[k] ;
                }
                if((neib[calc_id+15]<num_block && cy[k]==-1 && cz[k]==-1) || (neib[calc_id+18]<num_block && cy[k]==1 && cz[k]==1)){
                    cyy=cy[k] ; czz=cz[k] ;
                }
                if((neib[calc_id+16]<num_block && cy[k]==-1 && cz[k]==1) || (neib[calc_id+17]<num_block && cy[k]==1 && cz[k]==-1)){
                    cyy=cy[k] ; czz=cz[k] ;
                }

                if((neib[calc_id+11]>=num_block && cxx==-1 && czz==-1) || (neib[calc_id+14]>=num_block && cxx==1 && czz==1)){
                    cxx=0 ; czz=0 ;
                }
                if((neib[calc_id+12]>=num_block && cxx==-1 && czz==1) || (neib[calc_id+13]>=num_block && cxx==1 && czz==-1)){
                    cxx=0 ; czz=0 ;
                }
                if((neib[calc_id+15]>=num_block && cyy==-1 && czz==-1) || (neib[calc_id+18]>=num_block && cyy==1 && czz==1)){
                    cyy=0 ; czz=0 ;
                }
                if((neib[calc_id+16]>=num_block && cyy==-1 && czz==1) || (neib[calc_id+17]>=num_block && cyy==1 && czz==-1)){
                    cyy=0 ; czz=0 ;
                }
                }
            nextB.push_back( nei[abs(cxx)][abs(cyy)][abs(czz)][k]) ;
            nextK.push_back(kinv[abs(cxx)][abs(cyy)][abs(czz)][k]) ;
        }
    }
}

template<typename Typ>
__global__ void set_wall_propagation(Typ *items, int boudrary_condition, int *neib, int *nextK, int *nextB, Typ *ftmp){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k, k_inv[27]={0} ;  set_k_inv(k_inv,(int)items[IDX_Q]) ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        if(boudrary_condition==1){ // bounce back noslip
            for(k=1;k<items[IDX_Q];k++){
                ftmp[id_f+k] = ftmp[neib[id_f]*(int)items[IDX_Q]+k_inv[k]] ;
            }        
        }
        else if(boudrary_condition==0){ // slip
            int nextB_id ; 
            for(k=1;k<items[IDX_Q];k++){
                nextB_id = neib[id_f]*(int)items[IDX_Q] + nextB[(int)(id_rho-items[IDX_num_calc])*(int)items[IDX_Q]+k] ;
                ftmp[id_f+k] = ftmp[neib[nextB_id]*(int)items[IDX_Q] + nextK[(int)(id_rho-items[IDX_num_calc])*(int)items[IDX_Q]+k] ] ;
            }
        }
    }
}

template void set_walln<float>(const vector<float>& , const vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& ) ;
template void set_bound2D<float>(vector<float>& , int, const vector<int>& , vector<int>& , vector<int>& ) ;
template void set_bound3D<float>(vector<float>& , int, const vector<int>& , vector<int>& , vector<int>& ) ;
template __global__ void set_wall_propagation<float>(float* , int, int*, int*, int*, float*) ;

template void set_walln<double>(const vector<double>& , const vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& , vector<int>& ) ;
template void set_bound2D<double>(vector<double>& , int, const vector<int>& , vector<int>& , vector<int>& ) ;
template void set_bound3D<double>(vector<double>& , int, const vector<int>& , vector<int>& , vector<int>& ) ;
template __global__ void set_wall_propagation<double>(double* , int, int*, int *, int*, double*) ;
