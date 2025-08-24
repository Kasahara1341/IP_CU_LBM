#include "../all.hpp"

template<typename Typ>
void set_Lagrange_coeff(vector<Typ> &items_adv, double ratio){
    int i, j ;
    Typ posx[5] ; for(i=0;i<5;i++) posx[i] = i*ratio ;
    Typ x_alpha = 2.5 * ratio - 0.5 ;
    Typ tmp ;
    for(i=0;i<5;i++){
        tmp = 1.0 ;
        for(j=0;j<5;j++){
            if(i!=j){
                tmp *= (x_alpha - posx[j])/(posx[i] - posx[j]) ;
            }
        }
        items_adv[4-i] = tmp ;
    }
}

template<typename Typ>
__global__ void set_wall_f(Typ *items, int *neib, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho>=items[8] && id_rho<items[8]+items[9]){
        int k ; 
        for(k=0;k<items[6];k++){
            f[id_f+k] = f[neib[id_f]*(int)items[6]+k] ;
        }
    }
}


template<typename Typ>
__global__ void set_out(Typ *items, Typ *items_adv, int *neib, Typ *ftmp, Typ *fout, Typ *fin){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho<items[8]){
        int k, cp[27]={0}, cm[27]={0} ; set_cpm(cp,cm,(int)items[6]) ;
        int cxp[27]={0}, cxm[27]={0}, cyp[27]={0}, cym[27]={0} ; set_cxypm(cxp,cxm,cyp,cym) ;
        int id_cp, id_cpp, id_cm, id_cmm ;
        int num1, num2 ;
        int k_number1[27], k_number2[12] ;
        if((int)items[6]==9){
            num1=9 ; num2 = 0 ; for(k=0;k<9;k++){k_number1[k]=k;}
        }
        else{
            num1 = 15 ; num2 = 12 ;
            for(k=0;k<7;k++){
                k_number1[k]=k ; k_number1[k+7]= 11 + k ;
            }   k_number1[14]=18 ;
            for(k=0;k<4;k++){
                k_number2[k]=7+k ; k_number2[k+4]=19+k ; k_number2[k+8]=23+k ;
            }

            num1 = 27 ; num2=0 ; for(k=0;k<27;k++){k_number1[k]=k;}
        }
        for(k=0;k<num1;k++){ // monolinear
            id_cp   = neib[id_f  +cp[k_number1[k]]]*(int)items[6] ; id_cm   = neib[id_f  +cm[k_number1[k]]]*(int)items[6] ;
            id_cpp  = neib[id_cp +cp[k_number1[k]]]*(int)items[6] ; id_cmm  = neib[id_cm +cm[k_number1[k]]]*(int)items[6] ;
            fout[id_f+k_number1[k]] = (items_adv[0]*ftmp[id_cpp+k_number1[k]] + items_adv[1]*ftmp[id_cp+k_number1[k]] + items_adv[2]*ftmp[id_f+k_number1[k]] 
            + items_adv[3]*ftmp[id_cm+k_number1[k]] + items_adv[4]*ftmp[id_cmm+k_number1[k]])/items[7] ; 
            fin[id_f+k_number1[k]]  = fout[id_f+k_number1[k]] ; 
        }
        for(k=0;k<num2;k++){  // bilinear ?
            Typ flaxx, flaxy, flaxxy ;
            id_cp   = neib[id_f  +cxp[k_number2[k]]]*(int)items[6] ; id_cm   = neib[id_f  +cxm[k_number2[k]]]*(int)items[6] ;
            id_cpp  = neib[id_cp +cxp[k_number2[k]]]*(int)items[6] ; id_cmm  = neib[id_cm +cxm[k_number2[k]]]*(int)items[6] ;
            flaxx = (items_adv[0]*ftmp[id_cpp+k_number2[k]] + items_adv[1]*ftmp[id_cp+k_number2[k]] + items_adv[2]*ftmp[id_f+k_number2[k]] 
            + items_adv[3]*ftmp[id_cm+k_number2[k]] + items_adv[4]*ftmp[id_cmm+k_number2[k]])/items[7] ;

            id_cp   = neib[id_f  +cyp[k_number2[k]]]*(int)items[6] ; id_cm   = neib[id_f  +cym[k_number2[k]]]*(int)items[6] ;
            id_cpp  = neib[id_cp +cyp[k_number2[k]]]*(int)items[6] ; id_cmm  = neib[id_cm +cym[k_number2[k]]]*(int)items[6] ;
            flaxy = (items_adv[0]*ftmp[id_cpp+k_number2[k]] + items_adv[1]*ftmp[id_cp+k_number2[k]] + items_adv[2]*ftmp[id_f+k_number2[k]] 
                + items_adv[3]*ftmp[id_cm+k_number2[k]] + items_adv[4]*ftmp[id_cmm+k_number2[k]])/items[7] ;
    
            id_cp   = neib[id_f  +cp[k_number2[k]]]*(int)items[6] ; id_cm   = neib[id_f  +cm[k_number2[k]]]*(int)items[6] ;
            id_cpp  = neib[id_cp +cp[k_number2[k]]]*(int)items[6] ; id_cmm  = neib[id_cm +cm[k_number2[k]]]*(int)items[6] ;
            flaxxy = (items_adv[0]*ftmp[id_cpp+k_number2[k]] + items_adv[1]*ftmp[id_cp+k_number2[k]] + items_adv[2]*ftmp[id_f+k_number2[k]] 
                + items_adv[3]*ftmp[id_cm+k_number2[k]] + items_adv[4]*ftmp[id_cmm+k_number2[k]])/items[7] ;

            fout[id_f+k_number2[k]] = flaxx + flaxy - flaxxy/items[7] ; // (flaxx-flaxxy/ratio) + (flaxy-flaxxy/ratio) + flaxxy/ratio
            // fin[id_f+k_number2[k]] = fout[id_f+k_number2[k]] ;
            atomicAdd(&fin[id_f+k_number2[k]], flaxxy/items[7]) ;
            atomicAdd(&fin[neib[id_f  +cxm[k_number2[k]]]*(int)items[6]+k_number2[k]], flaxx-flaxxy/items[7]) ;
            atomicAdd(&fin[neib[id_f  +cym[k_number2[k]]]*(int)items[6]+k_number2[k]], flaxy-flaxxy/items[7]) ; // */
        } 
    }
}

template<typename Typ>
__global__ void set_wall_f_inv(Typ *items, int *neib, int *nextB, int *nextK, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho>=items[8] && id_rho<items[8]+items[9]){
        int k ; 
        int k_inv[27]={0} ; //set_k_inv(k_inv,(int)items[6]) ;
        if((int)items[6]==9){
            k_inv[1]=3 ; k_inv[3]=1 ; k_inv[5]=6 ; k_inv[6]=5 ; k_inv[7]=8 ; k_inv[8]=7 ;
            for(k=1;k<items[6];k++){
                f[id_f+k] = f[neib[id_f]*(int)items[6]+k_inv[k]] ;
            }
        }
        else if((int)items[6]==27){
            int k_number1[15], k_number2[12], num1, num2, nkz[8] ; 
            nkz[0 ]=20 ; nkz[1 ]=19 ; nkz[2 ]=22 ; nkz[3 ]=21 ;
            nkz[4 ]=24 ; nkz[5 ]=23 ; nkz[6 ]=26 ; nkz[7 ]=25 ;
            num1 = 15 ; num2 = 4 ;
            for(k=0;k<7;k++){
                k_number1[k]=k ; k_number1[k+7]= 11 + k ;
            }   k_number1[14]=18 ;
            for(k=0;k<4;k++){
                k_number2[k]=7+k ; k_number2[k+4]=19+k*2 ; k_number2[k+8]=20+k*2 ;
            }// 7 8 9 10 19 21 23 25 20 22 24 26
            set_k_inv(k_inv,(int)items[6]) ;
            int k_invxy[19] ;
            k_invxy[0 ]=0  ; k_invxy[1 ]=3  ; k_invxy[2 ]=4  ; k_invxy[3 ]=1  ; k_invxy[4 ]=2  ; k_invxy[5 ]=0  ; k_invxy[6 ]=0  ; k_invxy[7 ]=0  ; k_invxy[8 ]=0  ;
            k_invxy[9 ]=0  ; k_invxy[10]=0  ; k_invxy[11]=13 ; k_invxy[12]=14 ; k_invxy[13]=11 ; k_invxy[14]=12 ; k_invxy[15]=17 ; k_invxy[16]=18 ; k_invxy[17]=15 ; 
            k_invxy[18]=16 ; 

            for(k=0;k<num1;k++){
                int wall_id =  (id_rho - (int)items[8])*(int)items[6] ;
                int calc_id = neib[id_f]*(int)items[6] ;
                // int nk = k_inv[nextK[wall_id+k_inv[k_number1[k]]]], neB = k_invxy[nextB[wall_id+k_inv[k_number1[k]]]] ;
                // f[id_f+k_number1[k]] = f[neib[calc_id+neB]*(int)items[6]+nk] ;
                f[id_f+k_number1[k]] = f[calc_id+k_invxy[k_number1[k]]] ;
                if(id_rho==(int)items[8]+58+0){
                    // printf("k=%d  nextB=%d nextK=%d\n",k_number1[k],0,k_invxy[k_number1[k]]) ;
                }
            }
            for(k=0;k<num2;k++){
                int wall_id =  (id_rho - (int)items[8])*(int)items[6] ;
                int calc_id = neib[id_f]*(int)items[6] ;
                int nk = k_inv[nextK[wall_id+k_inv[k_number2[k]]]], neB = nextB[wall_id+k_inv[k_number2[k]]] ;
                f[id_f+k_number2[k+0]] = f[neib[calc_id+neB]*(int)items[6] + nk] ;
                f[id_f+k_number2[k+4]] = f[neib[calc_id+neB]*(int)items[6] + nkz[2*(nk-7) + 1]] ;
                f[id_f+k_number2[k+8]] = f[neib[calc_id+neB]*(int)items[6] + nkz[2*(nk-7) + 0]] ;
                if(id_rho==(int)items[8]){
                    // printf("k=%d  nextB=%d nextK=%d\n",k_number1[k],neB,nk) ;
                }
            }
        }
    }
}
template<typename Typ> // add fin from wall lattice 
__global__ void set_fin(Typ *items, int *neib, int *nextK, Typ *fin){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho>=items[8] && id_rho<items[8]+items[9]){
        if((int)items[6]==27){
            int k, l, k_inv[27]={0} ; set_k_inv(k_inv,(int)items[6]) ;
            int k_number2[12], nk[8], nkz[8] ; 
            nk[0]=9 ; nk[1]=8 ; nk[2]=10 ; nk[3]=7 ; nk[4]=7 ; nk[5]=10 ; nk[6]=8 ; nk[7]=9 ;
            nkz[0 ]=19 ; nkz[1 ]=20 ; nkz[2 ]=21 ; nkz[3 ]=22 ;
            nkz[4 ]=23 ; nkz[5 ]=24 ; nkz[6 ]=25 ; nkz[7 ]=26 ;
            int nneib[8] ; // 0,1 => 7 ; 2,3 => 8 ; 4,5 => 9 ; 6,7 => 10  x,y
            nneib[0]=1 ; nneib[1]=2 ; nneib[2]=1 ; nneib[3]=4 ; nneib[4]=3 ; nneib[5]=2 ; nneib[6]=3 ; nneib[7]=4 ;
            for(k=0;k<4;k++){
                k_number2[k]=7+k ; k_number2[k+4]=19+k*2 ; k_number2[k+8]=20+k*2 ;
            }// 7 8 9 10 19 21 23 25 20 22 24 26
            for(k=0;k<4;k++){ // devide by number of fin add to wall min=>0 max=>2
                fin[id_f + k_number2[k+0]] /= ((neib[id_f]+1)/(neib[id_f + nneib[k*2+0]]+1)) + ((neib[id_f]+1)/(neib[id_f + nneib[k*2+1]]+1)) + powf(10,-8) ; 
                fin[id_f + k_number2[k+4]] /= ((neib[id_f]+1)/(neib[id_f + nneib[k*2+0]]+1)) + ((neib[id_f]+1)/(neib[id_f + nneib[k*2+1]]+1)) + powf(10,-8) ; 
                fin[id_f + k_number2[k+8]] /= ((neib[id_f]+1)/(neib[id_f + nneib[k*2+0]]+1)) + ((neib[id_f]+1)/(neib[id_f + nneib[k*2+1]]+1)) + powf(10,-8) ; 
            }
            for(k=0;k<4;k++){ // add fin to neighbor lattice of wall
                for(l=0;l<2;l++){
                    fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + nk[ 2*(k_number2[k+0] -7) + l]]               += fin[id_f + k_number2[k+0]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;
                    fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + nkz[2*(nk[2*(k_number2[k+0] -7) + l]-7) + 0]] += fin[id_f + k_number2[k+4]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;
                    fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + nkz[2*(nk[2*(k_number2[k+0] -7) + l]-7) + 1]] += fin[id_f + k_number2[k+8]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;

                    // fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + k_number2[k+0]] += fin[id_f + k_number2[k+0]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;
                    // fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + k_number2[k+4]] += fin[id_f + k_number2[k+4]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;
                    // fin[neib[id_f+nneib[k*2+l]]*(int)items[6] + k_number2[k+8]] += fin[id_f + k_number2[k+8]] * (1-neib[id_f + nneib[k*2+l]]/id_rho) ;
                    if(id_rho==(int)items[8]+58*0+0){
                        // printf(" k_number2=%02d l=%d neib[id_f+nneib[k*2+l]]=%05d nk =%02d \n",k_number2[k+0],l,neib[id_f+nneib[k*2+l]],nk[2*(k_number2[k+0] -7) + l]) ;
                        // printf(" k_number2=%02d l=%d neib[id_f+nneib[k*2+l]]=%05d nzk=%02d \n",k_number2[k+4],l,neib[id_f+nneib[k*2+l]],nkz[2*(nk[2*(k_number2[k+0] -7) + l]-7) + 0]) ;
                        // printf(" k_number2=%02d l=%d neib[id_f+nneib[k*2+l]]=%05d nzk=%02d \n",k_number2[k+8],l,neib[id_f+nneib[k*2+l]],nkz[2*(nk[2*(k_number2[k+0] -7) + l]-7) + 1]) ;
                    }
                }
            }          
            // if(id_rho==(int)items[8]+58*0){printf("\n") ;  }
            /*for(k=0;k<4;k++){
                int wall_id =  (id_rho - (int)items[8])*(int)items[6] ;
                int calc_id = neib[id_f]*(int)items[6] ;
                int nk = k_inv[nextK[wall_id+k_inv[k_number2[k]]]];
                fin[calc_id+k_number2[k+0]] += fin[id_f + nk] ; 
                fin[calc_id+k_number2[k+4]] += fin[id_f + nkz[2*(nk-7) + 1 ]] ;
                fin[calc_id+k_number2[k+8]] += fin[id_f + nkz[2*(nk-7) + 0 ]] ;
                if(id_rho==(int)items[8]+58*0+0){
                    // printf("k_number2=%02d calc_id=%d neib[k_number2[k]]=%05d nk=%02d fin[nk]=%f \n",k_number2[k],neib[id_f],neib[calc_id+k_number2[k]],nk,fin[id_f+nk]) ;
                    // printf("k_number2=%02d nk=%02d fin[nk]=%f \n",k_number2[k+4],nkz[2*(nk-7)+1],fin[id_f + nkz[2*(nk-7) + 1 ]]) ;
                    // printf("k_number2=%02d nk=%02d fin[nk]=%f neib1=%d neib2=%d \n",k_number2[k+8],nkz[2*(nk-7)+0],fin[id_f + nkz[2*(nk-7) + 0 ]],neib[calc_id+6],neib[calc_id+3]) ;
                }
            }// */
        }
    }
}
template<typename Typ>
__global__ void set_tmp(Typ *items, int *neib, Typ *f, Typ *ftmp, Typ *fout, Typ *fin){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho<items[8]){
        int k, cp[27]={0}, cm[27]={0} ; set_cpm(cp,cm,(int)items[6]) ;
        for(k=0;k<items[6];k++){
                ftmp[id_f+k] = f[neib[id_f+cp[k]]*(int)items[6]+k] + (fin[id_f+k]-fout[neib[id_f+cp[k]]*(int)items[6]+k]) ;
        }
    }
}    

template void set_Lagrange_coeff(vector<float>& , double) ;
template __global__ void set_wall_f<float>(float*, int*, float*) ;
template __global__ void set_wall_f_inv<float>(float*, int*, int*, int*, float*) ;
template __global__ void set_out<float>(float*, float*, int*, float*, float*, float*) ;
template __global__ void set_fin<float>(float*, int*, int*, float*) ;
template __global__ void set_tmp<float>(float*, int*, float*, float*, float*, float*) ;

template void set_Lagrange_coeff(vector<double>& , double) ;
template __global__ void set_wall_f<double>(double*, int*, double*) ;
template __global__ void set_wall_f_inv<double>(double*, int*, int*, int*, double*) ;
template __global__ void set_out<double>(double*, double*, int*, double*, double*, double*) ;
template __global__ void set_fin<double>(double*, int*, int*, double*) ;
template __global__ void set_tmp<double>(double*, int*, double*, double*, double*, double*) ;
