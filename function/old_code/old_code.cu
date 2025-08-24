template<typename Typ>
__global__ void col_f_SRT(Typ *items, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    int k ; 
    if(id_rho<items[8]){
        for(k=0;k<items[6];k++){
            f[id_f+k] = f[id_f+k] - (f[id_f+k]-feq[id_f+k]+Fk[id_f+k]/2)/items[12] + Fk[id_f+k] ;
            ftmp[id_f+k] = f[id_f+k] ;
        }
    }
}
template<typename Typ>
__global__ void col_f_MRT(Typ *items, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk, Typ *M, Typ *M_inv, Typ *S){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    int k, l ; 
    if(id_rho<items[8]){
        for(k=0;k<items[6];k++){
            feq[id_f+k] = (f[id_f+k]-feq[id_f+k]+Fk[id_f+k]/2.0) ;
        }
        for(k=0;k<items[6];k++){
            Typ tmp=0 ;
            for(l=0;l<items[6];l++){
                tmp += M[k*(int)items[6]+l]*feq[id_f+l] ;
            }
            f[id_f+k]    = f[id_f+k] - tmp + Fk[id_f+k] ;
            ftmp[id_f+k] = f[id_f+k] ;
        }
    }
}

// QUICKESTが残っている
template<typename Typ>
__global__ void set_out(Typ *items, Typ *items_adv, int *neib, Typ *ftmp, Typ *fout){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho<items[8]){
        int k ; 
        int cp[27] ;
        int cm[27] ;
        if((int)items[6]==9){
            cp[0]=0 ; cp[1]=1 ; cp[2]=0 ; cp[3]=3 ; cp[4]=0 ; cp[5]=1 ; cp[6]=3 ; cp[7]=3 ; cp[8]=1 ;
            cm[0]=0 ; cm[1]=3 ; cm[2]=0 ; cm[3]=1 ; cm[4]=0 ; cm[5]=3 ; cm[6]=1 ; cm[7]=1 ; cm[8]=3 ;
        }else if((int)items[6]==27){
            cp[0 ]=0  ; cp[1 ]=1  ; cp[2 ]=2  ; cp[3 ]=3  ; cp[4 ]=4  ; cp[5 ]=0  ; cp[6 ]=0  ; cp[7 ]=7  ; cp[8 ]=8  ;
            cp[9 ]=9  ; cp[10]=10 ; cp[11]=1  ; cp[12]=1  ; cp[13]=3  ; cp[14]=3  ; cp[15]=2  ; cp[16]=2  ; cp[17]=4  ;
            cp[18]=4  ; cp[19]=7  ; cp[20]=7  ; cp[21]=8  ; cp[22]=8  ; cp[23]=9  ; cp[24]=9  ; cp[25]=10 ; cp[26]=10 ;
            cm[0 ]=0  ; cm[1 ]=3  ; cm[2 ]=4  ; cm[3 ]=1  ; cm[4 ]=2  ; cm[5 ]=0  ; cm[6 ]=0  ; cm[7 ]=10 ; cm[8 ]=9  ;
            cm[9 ]=8  ; cm[10]=7  ; cm[11]=3  ; cm[12]=3  ; cm[13]=1  ; cm[14]=1  ; cm[15]=4  ; cm[16]=4  ; cm[17]=2  ;
            cm[18]=2  ; cm[19]=10 ; cm[20]=10 ; cm[21]=9  ; cm[22]=9  ; cm[23]=8  ; cm[24]=8  ; cm[25]=7  ; cm[26]=7  ;
        }
        Typ coeff=1.0 ;//+ 0.008*items[7] ;
        for(k=0;k<items[6];k++){
            int id_cp = neib[id_f+cp[k]]*(int)items[6], id_cm = neib[id_f+cm[k]]*(int)items[6] ;
            int id_cpp = neib[id_cp+cp[k]]*(int)items[6], id_cmm = neib[id_cm+cm[k]]*(int)items[6] ;
            int id_cppp = neib[id_cpp+cp[k]]*(int)items[6], id_cmmm = neib[id_cmm+cm[k]]*(int)items[6] ;
            // fout[id_f+k] = ftmp[id_f+k] ; // 1st upwind }
            // QUICK
            // fout[id_f+k] = (3*ftmp[id_cp+k] + 6*ftmp[id_f+k] - ftmp[id_cm+k])/8.0 ;
            // Lagrange interpolation
            fout[id_f+k] = items_adv[0]*ftmp[id_cpp+k] + items_adv[1]*ftmp[id_cp+k] + items_adv[2]*ftmp[id_f+k] 
            + items_adv[3]*ftmp[id_cm+k] + items_adv[4]*ftmp[id_cmm+k] ; // */
            // QUICKEST TMU https://shintani.fpark.tmu.ac.jp/classes/open_university/advection_quickest.html
            //fout[id_f+k] = 0.5*(ftmp[id_f+k]+ftmp[id_cp+k]) - 0.5*(ftmp[id_cp+k]-ftmp[id_f+k])/items[7] - (1-1.0/pow(items[7],2))*(ftmp[id_cp+k]-2*ftmp[id_f+k]+ftmp[id_cm+k])/6.0 ; // QUICKEST
            // QUCIKEST https://www.jstage.jst.go.jp/article/kaigan/71/2/71_I_151/_pdf
            // fout[id_f+k] = (2*ftmp[id_cp+k] + 5*ftmp[id_f+k] - ftmp[id_cm+k])/6.0 - 0.5*(ftmp[id_cp+k]-ftmp[id_f+k])/items[7]/coeff + (ftmp[id_cp+k]-2*ftmp[id_f+k]+ftmp[id_cm+k])/(6.0*items[7]*items[7])/coeff ;
            // fout[id_f+k] = (3*ftmp[id_cp+k] + 6*ftmp[id_f+k] - ftmp[id_cm+k])/8.0 - 0.5*(ftmp[id_cp+k]-ftmp[id_f+k])/items[7]/coeff + (ftmp[id_cp+k]-2*ftmp[id_f+k]+ftmp[id_cm+k])/(6.0*items[7]*items[7])/coeff ;
            // fout[id_f+k] = (  ftmp[id_cp+k] +   ftmp[id_f+k]                )/2.0 - 0.5*(ftmp[id_cp+k]-ftmp[id_f+k])/items[7]/coeff + (ftmp[id_cp+k]-2*ftmp[id_f+k]+ftmp[id_cm+k])/(6.0*items[7]*items[7])/coeff ;
        }
    }
}


template<typename Typ>
__global__ void set_wall_f_inv(Typ *items, int *neib, int *nextB, int *nextK, Typ *f){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[6] ;
    if(id_rho>=items[8] && id_rho<items[8]+items[9]){
        int k ; 
        int z0[27], zeq[27], k_inv[27]={0} ; //set_k_inv(k_inv,(int)items[6]) ;
        if((int)items[6]==9){
            k_inv[1]=3 ; k_inv[3]=1 ; k_inv[5]=6 ; k_inv[6]=5 ; k_inv[7]=8 ; k_inv[8]=7 ;
            for(k=1;k<items[6];k++){
                f[id_f+k] = f[neib[id_f]*(int)items[6]+k_inv[k]] ;
            }
        }
        else if((int)items[6]==27){
            set_k_inv(k_inv,(int)items[6]) ;
            z0[0 ]=0  ; z0[1 ]=1  ; z0[2 ]=2  ; z0[3 ]=3  ; z0[4 ]=4  ; z0[5 ]=0  ; z0[6 ]=0  ; z0[7 ]=7  ; z0[8 ]=8  ;
            z0[9 ]=9  ; z0[10]=10 ; z0[11]=1  ; z0[12]=1  ; z0[13]=3  ; z0[14]=3  ; z0[15]=2  ; z0[16]=2  ; z0[17]=4  ; 
            z0[18]=4  ; z0[19]=7  ; z0[20]=7  ; z0[21]=8  ; z0[22]=8  ; z0[23]=9  ; z0[24]=9  ; z0[25]=10 ; z0[26]=10 ;
            zeq[0 ]=0  ; zeq[1 ]=1  ; zeq[2 ]=4  ; zeq[3 ]=3  ; zeq[4 ]=5  ; zeq[5 ]=0  ; zeq[6 ]=0  ; zeq[7 ]=7  ; zeq[8 ]=8  ;
            zeq[9 ]=9  ; zeq[10]=10 ; zeq[11]=11 ; zeq[12]=12 ; zeq[13]=13 ; zeq[14]=14 ; zeq[15]=15 ; zeq[16]=16 ; zeq[17]=17 ;
            zeq[18]=18 ; zeq[19]=19 ; zeq[20]=20 ; zeq[21]=21 ; zeq[22]=22 ; zeq[23]=23 ; zeq[24]=24 ; zeq[25]=25 ; zeq[26]=26 ;
            for(k=1;k<items[6];k++){
                int wall_id =  (id_rho - (int)items[8])*(int)items[6] ;
                int calc_id = neib[id_f]*(int)items[6] ;
                int nk = k_inv[nextK[wall_id+k_inv[k]]], neB = z0[nextB[wall_id+k_inv[k]]] ;
                f[id_f+k] = f[neib[calc_id+neB]*(int)items[6]+nk] ;
            }
        }
    }
}