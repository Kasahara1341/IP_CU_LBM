#include "../all.hpp"

template<typename Typ>
__global__ void set_wall_f(Typ *items, int *neib, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        int k ; 
        for(k=0;k<items[IDX_Q];k++){
            f[id_f+k] = f[neib[id_f]*(int)items[IDX_Q]+k] ;
        }
    }
}

template<typename Typ>
__device__ void set_posl_x(Typ *items, Typ *posl_x, int k){ // 中心から順番に壁への外挿か座標取得かを決定し計算する
    posl_x[1] = fminf(posl_x[1],posl_x[2]-items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[1],posl_x[2]-items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[0] = fminf(posl_x[0],posl_x[1]-items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[0],posl_x[1]-items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[3] = fmaxf(posl_x[3],posl_x[2]+items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fminf(posl_x[3],posl_x[2]+items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[4] = fmaxf(posl_x[4],posl_x[3]+items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fminf(posl_x[4],posl_x[3]+items[IDX_cx(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
}
template<typename Typ>
__device__ void set_posl_y(Typ *items, Typ *posl_x, int k){ // 中心から順番に壁への外挿か座標取得かを決定し計算する
    posl_x[1] = fminf(posl_x[1],posl_x[2]-items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[1],posl_x[2]-items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==-1.0) ;
    posl_x[0] = fminf(posl_x[0],posl_x[1]-items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[0],posl_x[1]-items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==-1.0) ;
    posl_x[3] = fmaxf(posl_x[3],posl_x[2]+items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==1.0) + fminf(posl_x[3],posl_x[2]+items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==-1.0) ;
    posl_x[4] = fmaxf(posl_x[4],posl_x[3]+items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==1.0) + fminf(posl_x[4],posl_x[3]+items[IDX_cy(k)]/items[IDX_c]*items[IDX_dz])*(items[IDX_cy(k)]/items[IDX_c]==-1.0) ;
}
template<typename Typ>
__device__ Typ distance(Typ posx0, Typ posy0, Typ posx1, Typ posy1){ // 中心から順番に壁への外挿か座標取得かを決定し計算する
    return sqrtf(powf(posx0-posx1,2)+powf(posy0-posy1,2)) ;
}
template<typename Typ>
__device__ void set_posl_xy(Typ *items, Typ *posl_x, int k){ // 中心から順番に壁への外挿か座標取得かを決定し計算する
    posl_x[1] = fminf(posl_x[1],posl_x[2]-items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[1],posl_x[2]-items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[0] = fminf(posl_x[0],posl_x[1]-items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fmaxf(posl_x[0],posl_x[1]-items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[3] = fmaxf(posl_x[3],posl_x[2]+items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fminf(posl_x[3],posl_x[2]+items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
    posl_x[4] = fmaxf(posl_x[4],posl_x[3]+items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==1.0) + fminf(posl_x[4],posl_x[3]+items[IDX_cx(k)]/items[IDX_c]*sqrtf(2)*items[IDX_dz])*(items[IDX_cx(k)]/items[IDX_c]==-1.0) ;
}
template<typename Typ>
__device__ void  Lagrange_interpolation(Typ *coeff, Typ *posl_x, Typ x_alpha){
    int i, j ; 
    for(i=0;i<5;i++){
        coeff[i] = 1.0 ;
        for(j=0  ;j<i;j++){coeff[i] *= (x_alpha - posl_x[j])/(posl_x[i]-posl_x[j]) ;}
        for(j=i+1;j<5;j++){coeff[i] *= (x_alpha - posl_x[j])/(posl_x[i]-posl_x[j]) ;}
    }
}

template<typename Typ> // test code
__global__ void set_out(Typ *items, int *neib, Typ *ftmp, Typ *fout, Typ *fin, Typ *finy, Typ *posx, Typ *posy, Typ *delX, Typ *delY){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        int k, cp[27]={0}, cm[27]={0} ; set_cpm(cp,cm,(int)items[IDX_Q]) ;
        int cxp[27]={0}, cxm[27]={0}, cyp[27]={0}, cym[27]={0} ; set_cxypm(cxp,cxm,cyp,cym) ;
        int id_cp, id_cpp, id_cm, id_cmm ;
        int num1=6, num2=6, num3=12 ;
        int k_number1[6], k_number2[6], k_number3[12] ;
        if((int)items[IDX_Q]==9){
            num2=0 ; num3 = 0 ;
            k_number1[0]=1 ; k_number1[1]=3 ; for(k=0;k<4;k++){k_number1[2+k]=5+k ;}
            fout[id_f]=ftmp[id_f] ; fout[id_f+2]=ftmp[id_f+2] ;fout[id_f+4]=ftmp[id_f+4] ;
            fin [id_f]=ftmp[id_f] ; fin [id_f+2]=ftmp[id_f+2] ;fin [id_f+4]=ftmp[id_f+4] ;
        }
        else{
            for(k=0;k<2;k++){
                k_number1[k]=1+k*2; k_number1[k+2]=11+k*2 ; k_number1[k+4]=12+k*2 ;
                k_number2[k]=2+k*2; k_number2[k+2]=15+k*2 ; k_number2[k+4]=16+k*2 ;
            }
            for(k=0;k<4;k++){ k_number3[k]=7+k ; k_number3[k+4]=19+k ; k_number3[k+8]=23+k ;}
            fout[id_f]=ftmp[id_f] ; fout[id_f+5]=ftmp[id_f+5] ; fout[id_f+6]=ftmp[id_f+6] ;
            fin [id_f]=ftmp[id_f] ; fin [id_f+5]=ftmp[id_f+5] ; fin [id_f+6]=ftmp[id_f+6] ;
        }
        for(k=0;k<num1;k++){ // monolinear x direction
            Typ coeff[5], posl_x[5], x_alpha = posx[id_rho] + items[IDX_cx(k_number1[k])]/items[IDX_c]*(delX[id_rho]-items[IDX_dz])/2.0 ;
            id_cp   = neib[id_f  +cp[k_number1[k]]]*(int)items[IDX_Q] ; id_cm   = neib[id_f  +cm[k_number1[k]]]*(int)items[IDX_Q] ;
            id_cpp  = neib[id_cp +cp[k_number1[k]]]*(int)items[IDX_Q] ; id_cmm  = neib[id_cm +cm[k_number1[k]]]*(int)items[IDX_Q] ;
            posl_x[0]=posx[id_cmm/(int)items[IDX_Q]] ; posl_x[1]=posx[id_cm/(int)items[IDX_Q]] ; posl_x[2]=posx[id_rho] ; posl_x[3]=posx[id_cp/(int)items[IDX_Q]] ; posl_x[4]=posx[id_cpp/(int)items[IDX_Q]] ;
            set_posl_x<Typ>(items,posl_x,k_number1[k]) ; // 計算点周りの座標決定が完了
            Lagrange_interpolation<Typ>(coeff,posl_x,x_alpha) ;
            fout[id_f+k_number1[k]] = (coeff[0]*ftmp[id_cmm+k_number1[k]] + coeff[1]*ftmp[id_cm+k_number1[k]] + coeff[2]*ftmp[id_f+k_number1[k]] 
            + coeff[3]*ftmp[id_cp+k_number1[k]] + coeff[4]*ftmp[id_cpp+k_number1[k]]) ;
            fin [id_f+k_number1[k]]  = fout[id_f+k_number1[k]]*(items[IDX_dz]/delX[int(id_cp/items[IDX_Q])]) ;
            fout[id_f+k_number1[k]] *= (items[IDX_dz]/delX[id_rho])  ;
        }
        for(k=0;k<num2;k++){ // monolinear y direction
            Typ coeff[5], posl_x[5], x_alpha = posy[id_rho] + items[IDX_cy(k_number2[k])]/items[IDX_c]*(delY[id_rho]-items[IDX_dz])/2.0 ;
            id_cp   = neib[id_f  +cp[k_number2[k]]]*(int)items[IDX_Q] ; id_cm   = neib[id_f  +cm[k_number2[k]]]*(int)items[IDX_Q] ;
            id_cpp  = neib[id_cp +cp[k_number2[k]]]*(int)items[IDX_Q] ; id_cmm  = neib[id_cm +cm[k_number2[k]]]*(int)items[IDX_Q] ;
            posl_x[0]=posy[id_cmm/(int)items[IDX_Q]] ; posl_x[1]=posy[id_cm/(int)items[IDX_Q]] ; posl_x[2]=posy[id_rho] ; posl_x[3]=posy[id_cp/(int)items[IDX_Q]] ; posl_x[4]=posy[id_cpp/(int)items[IDX_Q]] ;
            set_posl_y<Typ>(items,posl_x,k_number2[k]) ; // 計算点周りの座標決定が完了
            Lagrange_interpolation<Typ>(coeff,posl_x,x_alpha) ;
            fout[id_f+k_number2[k]] = (coeff[0]*ftmp[id_cmm+k_number2[k]] + coeff[1]*ftmp[id_cm+k_number2[k]] + coeff[2]*ftmp[id_f+k_number2[k]] 
            + coeff[3]*ftmp[id_cp+k_number2[k]] + coeff[4]*ftmp[id_cpp+k_number2[k]]) ;
            fin [id_f+k_number2[k]]  = fout[id_f+k_number2[k]]*(items[IDX_dz]/delY[int(id_cp/items[IDX_Q])]) ;
            fout[id_f+k_number2[k]] *=(items[IDX_dz]/delY[id_rho])  ;
        }
        for(k=0;k<num3;k++){  // bilinear ? 7,8,9,10 , 19~26
            Typ flaxx, flaxy, flaxxy ;
            int id_cpx, id_cpy, id_cpxy ;
            Typ coeff[5], posl_x[5], x_alpha ;
            x_alpha = posx[id_rho] + items[IDX_cx(k_number3[k])]/items[IDX_c]*(delX[id_rho]-items[IDX_dz])/2.0 ;
            id_cp   = neib[id_f  +cxp[k_number3[k]]]*(int)items[IDX_Q] ; id_cm   = neib[id_f  +cxm[k_number3[k]]]*(int)items[IDX_Q] ;
            id_cpp  = neib[id_cp +cxp[k_number3[k]]]*(int)items[IDX_Q] ; id_cmm  = neib[id_cm +cxm[k_number3[k]]]*(int)items[IDX_Q] ;
            posl_x[0]=posx[id_cmm/(int)items[IDX_Q]] ; posl_x[1]=posx[id_cm/(int)items[IDX_Q]] ; posl_x[2]=posx[id_rho] ; posl_x[3]=posx[id_cp/(int)items[IDX_Q]] ; posl_x[4]=posx[id_cpp/(int)items[IDX_Q]] ;
            set_posl_x<Typ>(items,posl_x,k_number3[k]) ; // 計算点周りの座標決定が完了
            Lagrange_interpolation<Typ>(coeff,posl_x,x_alpha) ;
            flaxx = (coeff[0]*ftmp[id_cmm+k_number3[k]] + coeff[1]*ftmp[id_cm+k_number3[k]] + coeff[2]*ftmp[id_f+k_number3[k]]
            + coeff[3]*ftmp[id_cp+k_number3[k]] + coeff[4]*ftmp[id_cpp+k_number3[k]]) ;
            id_cpx = id_cp/(int)items[IDX_Q] ;

            x_alpha = posy[id_rho] + items[IDX_cy(k_number3[k])]/items[IDX_c]*(delY[id_rho]-items[IDX_dz])/2.0 ;
            id_cp   = neib[id_f  +cyp[k_number3[k]]]*(int)items[IDX_Q] ; id_cm   = neib[id_f  +cym[k_number3[k]]]*(int)items[IDX_Q] ;
            id_cpp  = neib[id_cp +cyp[k_number3[k]]]*(int)items[IDX_Q] ; id_cmm  = neib[id_cm +cym[k_number3[k]]]*(int)items[IDX_Q] ;
            posl_x[0]=posy[id_cmm/(int)items[IDX_Q]] ; posl_x[1]=posy[id_cm/(int)items[IDX_Q]] ; posl_x[2]=posy[id_rho] ; posl_x[3]=posy[id_cp/(int)items[IDX_Q]] ; posl_x[4]=posy[id_cpp/(int)items[IDX_Q]] ;
            set_posl_y<Typ>(items,posl_x,k_number3[k]) ; // 計算点周りの座標決定が完了
            Lagrange_interpolation<Typ>(coeff,posl_x,x_alpha) ;
            flaxy = (coeff[0]*ftmp[id_cmm+k_number3[k]] + coeff[1]*ftmp[id_cm+k_number3[k]] + coeff[2]*ftmp[id_f+k_number3[k]] 
            + coeff[3]*ftmp[id_cp+k_number3[k]] + coeff[4]*ftmp[id_cpp+k_number3[k]]) ;
            id_cpy = id_cp/(int)items[IDX_Q] ;            

            id_cp   = neib[id_f  +cp[k_number3[k]]]*(int)items[IDX_Q] ; id_cm   = neib[id_f  +cm[k_number3[k]]]*(int)items[IDX_Q] ;
            id_cpp  = neib[id_cp +cp[k_number3[k]]]*(int)items[IDX_Q] ; id_cmm  = neib[id_cm +cm[k_number3[k]]]*(int)items[IDX_Q] ;
            posl_x[2] = 0 ; x_alpha = items[IDX_cx(k_number3[k])]/items[IDX_c]*(sqrtf(powf(delX[id_rho],2)+powf(delY[id_rho],2))-sqrtf(2.0)*items[IDX_dz])/2.0 ;
            posl_x[0] = -distance(posx[id_rho],posy[id_rho],posx[id_cmm/(int)items[IDX_Q]],posy[id_cmm/(int)items[IDX_Q]]) * items[IDX_cx(k_number3[k])]/items[IDX_c] ;
            posl_x[1] = -distance(posx[id_rho],posy[id_rho],posx[id_cm /(int)items[IDX_Q]],posy[id_cm /(int)items[IDX_Q]]) * items[IDX_cx(k_number3[k])]/items[IDX_c] ;
            posl_x[3] =  distance(posx[id_rho],posy[id_rho],posx[id_cp /(int)items[IDX_Q]],posy[id_cp /(int)items[IDX_Q]]) * items[IDX_cx(k_number3[k])]/items[IDX_c] ;
            posl_x[4] =  distance(posx[id_rho],posy[id_rho],posx[id_cpp/(int)items[IDX_Q]],posy[id_cpp/(int)items[IDX_Q]]) * items[IDX_cx(k_number3[k])]/items[IDX_c] ;
            set_posl_xy<Typ>(items,posl_x,k_number3[k]) ; // 計算点周りの座標決定が完了
            Lagrange_interpolation<Typ>(coeff,posl_x,x_alpha) ;
            flaxxy = (coeff[0]*ftmp[id_cmm+k_number3[k]] + coeff[1]*ftmp[id_cm+k_number3[k]] + coeff[2]*ftmp[id_f+k_number3[k]] 
            + coeff[3]*ftmp[id_cp+k_number3[k]] + coeff[4]*ftmp[id_cpp+k_number3[k]]) ;
            id_cpxy = id_cp/(int)items[IDX_Q] ;
            // fout[id_f+k_number3[k]] = flaxx*items[IDX_dz]/delX[id_rho] + flaxy*items[IDX_dz]/delY[id_rho] - flaxxy*(items[IDX_dz]*items[IDX_dz])/(delX[id_rho]*delY[id_rho]) ; 

            fout[id_f+k_number3[k]] = (flaxx - flaxxy*items[IDX_dz]/delY[id_rho])*items[IDX_dz]/delX[id_rho] * (1-int(items[IDX_dz]/delY[id_rho])) 
                                    + (flaxy - flaxxy*items[IDX_dz]/delX[id_rho])*items[IDX_dz]/delY[id_rho] * (1-int(items[IDX_dz]/delX[id_rho])) 
                                    +  flaxxy*(items[IDX_dz]*items[IDX_dz])/(delX[id_rho]*delY[id_rho]) ; // */ 

            atomicAdd(&fin[id_f+k_number3[k]], flaxxy*(items[IDX_dz]*items[IDX_dz])/(delX[id_cpx]*delY[id_cpy])) ;
            atomicAdd(&fin[neib[id_f   +cxm[k_number3[k]]]*(int)items[IDX_Q]+k_number3[k]], (flaxy-flaxxy*items[IDX_dz]/delX[id_cpy])*items[IDX_dz]/delY[id_cpy] * Typ(1-int(items[IDX_dz]/delX[id_rho])) ) ;
            atomicAdd(&finy[neib[id_f  +cym[k_number3[k]]]*(int)items[IDX_Q]+k_number3[k]], (flaxx-flaxxy*items[IDX_dz]/delY[id_cpx])*items[IDX_dz]/delX[id_cpx] * Typ(1-int(items[IDX_dz]/delY[id_rho])) ) ; // */

            // atomicAdd(&fin[neib[id_f   +cxm[k_number3[k]]]*(int)items[IDX_Q]+k_number3[k]], (flaxy-flaxy*items[IDX_dz]/delX[id_cpy])*items[IDX_dz]/delY[id_cpy]);// * (1-int(items[IDX_dz]/delX[id_rho])) ) ;
            // atomicAdd(&finy[neib[id_f  +cym[k_number3[k]]]*(int)items[IDX_Q]+k_number3[k]], (flaxx-flaxx*items[IDX_dz]/delY[id_cpx])*items[IDX_dz]/delX[id_cpx]);// * (1-int(items[IDX_dz]/delY[id_rho])) ) ; // */
            /*if(1-int(items[IDX_dz]/delX[id_rho]) == 1){
                printf("id=rho=%d dz=%e, delX=%e, atomicAdd=%f\n",id_rho,items[IDX_dz],delX[id_rho], 
                float(1-int(items[IDX_dz]/delX[id_rho])));
            } // */
        } 
    }
}
template<typename Typ> // add fin from wall lattice 
__global__ void set_fin(Typ *items, int *neib, int *nextK, Typ *fin, Typ *finy){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        if((int)items[IDX_Q]==27){
            int k, k_inv[27]={0} ; set_k_inv(k_inv,(int)items[IDX_Q]) ;
            int k_number3[12], nk[8], nkz[8] ; 
            nk[0]=9 ; nk[1]=8 ; nk[2]=10 ; nk[3]=7 ; nk[4]=7 ; nk[5]=10 ; nk[6]=8 ; nk[7]=9 ;
            nkz[0 ]=19 ; nkz[1 ]=20 ; nkz[2 ]=21 ; nkz[3 ]=22 ;
            nkz[4 ]=23 ; nkz[5 ]=24 ; nkz[6 ]=25 ; nkz[7 ]=26 ;
            int nneib[8] ; // 0,1 => 7 ; 2,3 => 8 ; 4,5 => 9 ; 6,7 => 10  x,y
            nneib[0]=1 ; nneib[1]=2 ; nneib[2]=1 ; nneib[3]=4 ; nneib[4]=3 ; nneib[5]=2 ; nneib[6]=3 ; nneib[7]=4 ;
            for(k=0;k<4;k++){
                k_number3[k]=7+k ; k_number3[k+4]=19+k*2 ; k_number3[k+8]=20+k*2 ;
            }// 7 8 9 10 19 21 23 25 20 22 24 26
            for(k=0;k<4;k++){ // add fin to neighbor lattice of wall
                fin[neib[id_f+nneib[k*2+0]]*(int)items[IDX_Q] + nk[ 2*(k_number3[k+0] -7) + 0]]               += fin[id_f + k_number3[k+0]] * (1-neib[id_f + nneib[k*2+0]]/id_rho) ;
                fin[neib[id_f+nneib[k*2+0]]*(int)items[IDX_Q] + nkz[2*(nk[2*(k_number3[k+0] -7) + 0]-7) + 0]] += fin[id_f + k_number3[k+4]] * (1-neib[id_f + nneib[k*2+0]]/id_rho) ;
                fin[neib[id_f+nneib[k*2+0]]*(int)items[IDX_Q] + nkz[2*(nk[2*(k_number3[k+0] -7) + 0]-7) + 1]] += fin[id_f + k_number3[k+8]] * (1-neib[id_f + nneib[k*2+0]]/id_rho) ;

                finy[neib[id_f+nneib[k*2+1]]*(int)items[IDX_Q] + nk[ 2*(k_number3[k+0] -7) + 1]]               += finy[id_f + k_number3[k+0]] * (1-neib[id_f + nneib[k*2+1]]/id_rho) ;
                finy[neib[id_f+nneib[k*2+1]]*(int)items[IDX_Q] + nkz[2*(nk[2*(k_number3[k+0] -7) + 1]-7) + 0]] += finy[id_f + k_number3[k+4]] * (1-neib[id_f + nneib[k*2+1]]/id_rho) ;
                finy[neib[id_f+nneib[k*2+1]]*(int)items[IDX_Q] + nkz[2*(nk[2*(k_number3[k+0] -7) + 1]-7) + 1]] += finy[id_f + k_number3[k+8]] * (1-neib[id_f + nneib[k*2+1]]/id_rho) ;
            }          
        }
    }
}
template<typename Typ>
__global__ void reset_f(Typ *items, Typ *ftmp, Typ *fin, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        int k ;
        for(k=0;k<items[IDX_Q];k++){
            fin[id_f+k] += f[id_f+k] ;
            f[id_f+k]   =  ftmp[id_f+k] ;
        }
    }
}
template<typename Typ>
__global__ void set_wall_f_inv(Typ *items, int *neib, int *nextB, int *nextK, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        int k ; 
        int k_inv[27]={0} ; //set_k_inv(k_inv,(int)items[IDX_Q]) ;
        if((int)items[IDX_Q]==9){
            k_inv[1]=3 ; k_inv[3]=1 ; k_inv[5]=6 ; k_inv[6]=5 ; k_inv[7]=8 ; k_inv[8]=7 ;
            for(k=1;k<items[IDX_Q];k++){
                f[id_f+k] = f[neib[id_f]*(int)items[IDX_Q]+k_inv[k]] ;
            }
        }
        else if((int)items[IDX_Q]==27){
            int k_number1[15], k_number3[12], num1, num3, nkz[8] ; 
            nkz[0 ]=20 ; nkz[1 ]=19 ; nkz[2 ]=22 ; nkz[3 ]=21 ;
            nkz[4 ]=24 ; nkz[5 ]=23 ; nkz[6 ]=26 ; nkz[7 ]=25 ;
            num1 = 15 ; num3 = 4 ;
            for(k=0;k<7;k++){
                k_number1[k]=k ; k_number1[k+7]= 11 + k ;
            }   k_number1[14]=18 ;
            for(k=0;k<4;k++){
                k_number3[k]=7+k ; k_number3[k+4]=19+k*2 ; k_number3[k+8]=20+k*2 ;
            }// 7 8 9 10 19 21 23 25 20 22 24 26
            set_k_inv(k_inv,(int)items[IDX_Q]) ;
            int k_invxy[19] ;
            k_invxy[0 ]=0  ; k_invxy[1 ]=3  ; k_invxy[2 ]=4  ; k_invxy[3 ]=1  ; k_invxy[4 ]=2  ; k_invxy[5 ]=0  ; k_invxy[6 ]=0  ; k_invxy[7 ]=0  ; k_invxy[8 ]=0  ;
            k_invxy[9 ]=0  ; k_invxy[10]=0  ; k_invxy[11]=13 ; k_invxy[12]=14 ; k_invxy[13]=11 ; k_invxy[14]=12 ; k_invxy[15]=17 ; k_invxy[16]=18 ; k_invxy[17]=15 ; 
            k_invxy[18]=16 ; 

            for(k=0;k<num1;k++){
                int calc_id = neib[id_f]*(int)items[IDX_Q] ;
                f[id_f+k_number1[k]] = f[calc_id+k_invxy[k_number1[k]]] ;
            }
            for(k=0;k<num3;k++){
                int wall_id =  (id_rho - (int)items[IDX_num_calc])*(int)items[IDX_Q] ;
                int calc_id = neib[id_f]*(int)items[IDX_Q] ;
                int nk = k_inv[nextK[wall_id+k_inv[k_number3[k]]]], neB = nextB[wall_id+k_inv[k_number3[k]]] ;
                f[id_f+k_number3[k+0]] = f[neib[calc_id+neB]*(int)items[IDX_Q] + nk] ;
                f[id_f+k_number3[k+4]] = f[neib[calc_id+neB]*(int)items[IDX_Q] + nkz[2*(nk-7) + 1]] ;
                f[id_f+k_number3[k+8]] = f[neib[calc_id+neB]*(int)items[IDX_Q] + nkz[2*(nk-7) + 0]] ;
            }
        }
    }
}

template<typename Typ>
__global__ void set_tmp(Typ *items, int *neib, Typ *f, Typ *ftmp, Typ *fout, Typ *fin){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        int k, cp[27]={0}, cm[27]={0} ; set_cpm(cp,cm,(int)items[IDX_Q]) ;
        for(k=0;k<items[IDX_Q];k++){
                ftmp[id_f+k] = f[neib[id_f+cp[k]]*(int)items[IDX_Q]+k] + (fin[id_f+k]-fout[neib[id_f+cp[k]]*(int)items[IDX_Q]+k]) ;
        }
    }
}
template<typename Typ>
void IP_process(Typ *d_items, int numBlocks, int blockSize, int *d_neib, Typ *d_f, Typ *d_feq, Typ *d_ftmp, Typ *d_fout, int *d_nextB, int *d_nextK, Typ *d_posx, Typ *d_posy, Typ* d_delX, Typ *d_delY, int slip){        
        set_wall_f    <Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_ftmp) ;
        set_out       <Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_ftmp, d_fout, d_feq, d_f, d_posx, d_posy, d_delX, d_delY) ;
        set_fin       <Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextK, d_feq, d_f) ;
        reset_f       <Typ> <<<numBlocks, blockSize>>>(d_items, d_ftmp, d_feq, d_f) ;
        set_wall_f_inv<Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_f) ;
        set_wall_f_inv<Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_fout) ;
        set_tmp       <Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_f, d_ftmp, d_fout, d_feq) ; // */
        set_wall_propagation<Typ> <<<numBlocks, blockSize>>>(d_items, slip, d_neib, d_nextK, d_nextB, d_ftmp) ; 
        propagation         <Typ> <<<numBlocks, blockSize>>>(d_items, d_neib, d_f, d_ftmp) ;
        reset_wall_fin      <Typ> <<<numBlocks, blockSize>>>(d_items, d_feq, d_f) ; 
}

template void IP_process<float>(float*, int, int, int*, float*, float*, float*, float*, int*, int*, float*, float*, float*, float*, int) ;
template void IP_process<double>(double*, int, int, int*, double*, double*, double*, double*, int*, int*, double*, double*, double*, double*, int) ;
