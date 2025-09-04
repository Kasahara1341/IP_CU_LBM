#include "../all.hpp"

template<typename Typ>
__global__ void equ_f(Typ *items, Typ *feq, Typ *pressure, Typ *velx, Typ *vely, Typ *velz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; Typ tmp ;
    if(id_rho<items[IDX_num_calc]){
        for(k=0;k<items[IDX_Q];k++){
            tmp = (items[IDX_cx(k)]*velx[id_rho]  + items[IDX_cy(k)]*vely[id_rho] +  items[IDX_cz(k)]*velz[id_rho])/(items[IDX_c]*items[IDX_c]) ;
            feq [id_f+k]  = items[IDX_w(k)] * (pressure[id_rho] + 3.0*tmp + 4.5*tmp*tmp-1.5*(velx[id_rho]*velx[id_rho]+vely[id_rho]*vely[id_rho]+velz[id_rho]*velz[id_rho])/(items[IDX_c]*items[IDX_c])) ;
        }
    }
}

template<typename Typ>
__global__ void equ_g(Typ *items, Typ *feq, Typ *sal, Typ *velx, Typ *vely, Typ *velz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; Typ tmp ;
    if(id_rho<items[IDX_num_calc]){
        for(k=0;k<items[IDX_Q];k++){
            tmp = (items[IDX_cx(k)]*velx[id_rho]  + items[IDX_cy(k)]*vely[id_rho] +  items[IDX_cz(k)]*velz[id_rho])/(items[IDX_c]*items[IDX_c]) ;
            feq [id_f+k]  = items[IDX_w(k)]*sal[id_rho] * (1.0 + 3.0*tmp + 4.5*tmp*tmp-1.5*(velx[id_rho]*velx[id_rho]+vely[id_rho]*vely[id_rho]+velz[id_rho]*velz[id_rho])/(items[IDX_c]*items[IDX_c])) ;
        }
    }
}
template<typename Typ>
__global__ void wall_function(Typ *items, Typ *delX, Typ *delY, int x, int y, int z, int wall_number, int *wall, Typ *velx, Typ *vely, Typ *velz, Typ *Fx, Typ *Fy, Typ *rho){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<wall_number ){ 
        int id=wall[id_rho] ;
        Typ volume = delX[id]*delY[id]*items[IDX_dz] ;
        Typ abs_u = sqrtf(velx[id]*velx[id] + vely[id]*vely[id] + velz[id]*velz[id]) ;
        Typ eta= abs_u*powf(10,-7) + powf(10,-20) ;
        Typ utau, epsilon, kap=0.41, Beta=5.5, tmp ;
        Typ y_dist = x*delX[id]/2.0 + y*delY[id]/2.0 + z*items[IDX_dz]/2.0 ;
        Typ area = volume/(y_dist*2.0) ; 
        utau=abs_u*0.11 ;
        // utau= sqrtf(abs_u*items[IDX_nu]/y_dist) ;
        for(int i=0;i<15;i++){
            tmp = kap*abs_u/(utau+eta) ;
            // tmp = fminf(8,tmp) ; // exp(tmp)のオーバーフローを防ぐため最大値を10とする
            /* calculate changing value */
            epsilon=
            -(
                (abs_u/(utau+eta) - y_dist*utau/items[IDX_nu] + exp(-kap*Beta)*( exp(tmp) - 1.0 - tmp - powf(tmp,2)/2.0 - powf(tmp,3)/6.0 )) /  // function of spalding
                (-y_dist/items[IDX_nu] - abs_u/(powf(utau,2)+eta) + exp(-kap*Beta)*tmp/(utau+eta)*( -exp(tmp) + 1.0 + tmp + powf(tmp,2)/(2.0) ) + eta)  // first derivative of above function
            ) ;
            utau+=epsilon ;    
        }
        Fx[id] += -velx[id]/(abs_u+eta)*powf(utau,2.0)* rho[id] * (area/volume) ;
        Fy[id] += -vely[id]/(abs_u+eta)*powf(utau,2.0)* rho[id] * (area/volume) ;
    }
}

template<typename Typ>
__global__ void Force(Typ *items, bool Boussi_flag, int *neib, Typ *f, Typ *feq, Typ *tau, Typ *Fk, Typ *Fx, Typ *Fy, Typ *Fz, Typ *pressure, Typ *rho, Typ *sal, Typ *phi, Typ *velx, Typ *vely, Typ *velz, Typ *delX, Typ *delY, Typ *posx, Typ *posy, Typ *posz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; 
    if(id_rho<items[IDX_num_calc]){
        int up=2, down=4 ;
        if((int)items[IDX_Q]==27){ up=5; down=6;}
        Typ nablax = 0, nablay = 0, nablaz = 0, nabla2 = 0 ;
        nablax = (phi[neib[id_f+1]]-phi[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        if((int)items[IDX_Q]==27){
        nablay = (phi[neib[id_f+2]]-phi[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;}
        nablaz = (phi[neib[id_f+up]]-phi[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;
        nabla2 = (phi[neib[id_f+1]] - 2*phi[id_rho] + phi[neib[id_f+3]])/powf(delX[id_rho],2) 
               + (phi[neib[id_f+2]] - 2*phi[id_rho] + phi[neib[id_f+4]])/powf(items[IDX_dz],2) ;
        // surface tension
        Typ beta, kap ; beta = 12.0*items[IDX_sigma]/items[IDX_PFthick] ; kap = 3.0*items[IDX_sigma]*items[IDX_PFthick]/2.0 ;
        Fx[id_rho] += (4*beta*phi[id_rho]*(phi[id_rho]-1.0)*(phi[id_rho]-0.5) - kap*nabla2) * nablax ;
        Fy[id_rho] += (4*beta*phi[id_rho]*(phi[id_rho]-1.0)*(phi[id_rho]-0.5) - kap*nabla2) * nablay ;
        Fz[id_rho] += (4*beta*phi[id_rho]*(phi[id_rho]-1.0)*(phi[id_rho]-0.5) - kap*nabla2) * nablaz ; // */
        // body force
        if(Boussi_flag){Fz[id_rho]+= - sal[id_rho] * 0.824493 * 9.81 ;}
        else{Fz[id_rho] += - rho[id_rho] * 9.81 ;}
        // pressure gradient
        nablax = 0 ; nablay = 0 ; nablaz = 0 ;
        nablax = (rho[neib[id_f+1]]-rho[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        if((int)items[IDX_Q]==27){
        nablay = (rho[neib[id_f+2]]-rho[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;}
        nablaz = (rho[neib[id_f+up]]-rho[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;
        Fx[id_rho] += - pressure[id_rho] * items[IDX_c]*items[IDX_c]/3.0 * nablax ;
        Fy[id_rho] += - pressure[id_rho] * items[IDX_c]*items[IDX_c]/3.0 * nablay ;
        Fz[id_rho] += - pressure[id_rho] * items[IDX_c]*items[IDX_c]/3.0 * nablaz ;

        /*for(k=1;k<items[IDX_Q];k++){
            // Fmu
            Fx[id_rho]+=-(1.0-0.5/tau[id_rho])*items[IDX_cx(k)]*items[IDX_cy(k)]*(f[id_f+k]-feq[id_f+k])*nablay
                        -(1.0-0.5/tau[id_rho])*items[IDX_cx(k)]*items[IDX_cz(k)]*(f[id_f+k]-feq[id_f+k])*nablaz 
                        -(1.0-0.5/tau[id_rho])*items[IDX_cx(k)]*items[IDX_cx(k)]*(f[id_f+k]-feq[id_f+k])*nablax ;
            Fy[id_rho]+=-(1.0-0.5/tau[id_rho])*items[IDX_cy(k)]*items[IDX_cx(k)]*(f[id_f+k]-feq[id_f+k])*nablax
                        -(1.0-0.5/tau[id_rho])*items[IDX_cy(k)]*items[IDX_cz(k)]*(f[id_f+k]-feq[id_f+k])*nablaz 
                        -(1.0-0.5/tau[id_rho])*items[IDX_cy(k)]*items[IDX_cy(k)]*(f[id_f+k]-feq[id_f+k])*nablay ;
            Fz[id_rho]+=-(1.0-0.5/tau[id_rho])*items[IDX_cz(k)]*items[IDX_cx(k)]*(f[id_f+k]-feq[id_f+k])*nablax
                        -(1.0-0.5/tau[id_rho])*items[IDX_cz(k)]*items[IDX_cy(k)]*(f[id_f+k]-feq[id_f+k])*nablay 
                        -(1.0-0.5/tau[id_rho])*items[IDX_cz(k)]*items[IDX_cz(k)]*(f[id_f+k]-feq[id_f+k])*nablaz ;
        } // */
        for(k=0;k<items[IDX_Q];k++){
            Typ tmp = 3.0*(items[IDX_cx(k)]*velx[id_rho] + items[IDX_cy(k)]*vely[id_rho] + items[IDX_cz(k)]*velz[id_rho] )/(items[IDX_c]*items[IDX_c]) ;
            Fk[id_f+k] = items[IDX_w(k)] * items[IDX_dt] * 
            (items[IDX_cx(k)]*Fx[id_rho] + items[IDX_cy(k)]*Fy[id_rho] + items[IDX_cz(k)]*Fz[id_rho])
            // ( (items[IDX_cx(k)]-velx[id_rho] + tmp*items[IDX_cx(k)])*Fx[id_rho]
            // + (items[IDX_cy(k)]-vely[id_rho] + tmp*items[IDX_cy(k)])*Fy[id_rho]
            // + (items[IDX_cz(k)]-velz[id_rho] + tmp*items[IDX_cz(k)])*Fz[id_rho])
            / (rho[id_rho]*items[IDX_c]*items[IDX_c]/3.0) ; 
        }
    }
}


template<typename Typ>
__global__ void col_f_MRT(Typ *items, Typ *tau, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk, Typ *M, Typ *M_inv, Typ *S, Typ *Mconst){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int i, k, l ; 
    if(id_rho<items[IDX_num_calc]){
        Typ M_LES[27] ;
        for(k=0;k<items[IDX_Q];k++){
            feq[id_f+k] = (f[id_f+k]-feq[id_f+k]+Fk[id_f+k]/2.0) ;
        }
        for(i=0;i<items[IDX_Q];i++){
            Typ tmp=0 ;
            for(k=0;k<items[IDX_Q];k++){
                M_LES[k] = 0 ;
                if((int)items[IDX_Q]==27){
                    for(l=0;l<5;l++){
                        M_LES[k] += M_inv[i*(int)items[IDX_Q]+(l+5)] * M[(l+5)*(int)items[IDX_Q] + k] / tau[id_rho] ;
                    }
                }
                else if((int)items[IDX_Q]==9){
                    for(l=0;l<2;l++){
                        M_LES[k] += M_inv[i*(int)items[IDX_Q]+(l+7)] * M[(l+7)*(int)items[IDX_Q] + k] / tau[id_rho] ;
                    }
                }
            }            
            for(k=0;k<items[IDX_Q];k++){
                tmp += (Mconst[i*(int)items[IDX_Q]+k]+M_LES[k])*feq[id_f+k] ;
            }
            f[id_f+i]    = f[id_f+i] - tmp + Fk[id_f+i] ;
            ftmp[id_f+i] = f[id_f+i] ; 
        }
        for(k=0;k<items[IDX_Q];k++){
            feq[id_f+k] = 0 ; f[id_f+k] = 0 ;
        }
    }
}
template<typename Typ>
__global__ void col_f_SRT(Typ *items, Typ *tau, Typ *f, Typ *ftmp, Typ *feq, Typ *Fk){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k; 
    if(id_rho<items[IDX_num_calc]){
        for(k=0;k<items[IDX_Q];k++){
            f[id_f+k] = f[id_f+k] - (f[id_f+k]-feq[id_f+k]+Fk[id_f+k]/2.0)/tau[id_rho] + Fk[id_f+k] ;
            ftmp [id_f+k] = f[id_f+k] ;
        }
        for(k=0;k<items[IDX_Q];k++){
            feq[id_f+k] = 0 ; f[id_f+k] = 0 ;
        }
    }
}

template<typename Typ>
__global__ void col_g_reg(Typ *items, Typ *taus, Typ *g, Typ *gtmp, Typ *geq, Typ *sal, Typ *velx, Typ *vely, Typ *velz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; 
    if(id_rho<items[IDX_num_calc]){
        Typ gkx=0, gky=0, gkz=0 ;
        for(k=0;k<items[IDX_Q];k++){
            float tmp = (items[IDX_cx(k)]*velx[id_rho]  + items[IDX_cy(k)]*vely[id_rho] +  items[IDX_cz(k)]*velz[id_rho])/(items[IDX_c]*items[IDX_c]) ;
            geq [id_f+k]  = items[IDX_w(k)]*sal[id_rho] * (1.0 + 3.0*tmp + 4.5*tmp*tmp-1.5*(velx[id_rho]*velx[id_rho]+vely[id_rho]*vely[id_rho]+velz[id_rho]*velz[id_rho])/(items[IDX_c]*items[IDX_c])) ;
            gtmp[id_f+k] = g[id_f+k] - geq[id_f+k] ;
            gkx+= gtmp[id_f+k]*items[IDX_cx(k)] ; gky+= gtmp[id_f+k]*items[IDX_cy(k)] ; gkz+= gtmp[id_f+k]*items[IDX_cz(k)] ;
        }
        for(k=0;k<items[IDX_Q];k++){
            gtmp[id_f+k] = items[IDX_w(k)]*3
              *(gkx*items[IDX_cx(k)] + gky*items[IDX_cy(k)] + gkz*items[IDX_cz(k)])
              /(items[IDX_c]*items[IDX_c]) ;
            g[id_f+k] = geq[id_f+k] + (1.0-1.0/taus[id_rho])*gtmp[id_f+k] ;
            gtmp[id_f+k] = g[id_f+k] ; geq[id_f+k] = 0 ; g[id_f+k] = 0 ;
        } 
    }
}

template<typename Typ>
__global__ void col_PF(Typ *items, int *neib, Typ *taus, Typ *g, Typ *gtmp, Typ *geq, Typ *phi, Typ *velx, Typ *vely, Typ *velz, Typ *phiold, Typ *uold, Typ *vold, Typ *wold, Typ *posx, Typ *posy, Typ *posz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; 
    if(id_rho<items[IDX_num_calc]){
        int up=2, down=4 ; if((int)items[IDX_Q]==27){ up=5; down=6;}
        Typ nablax = 0, nablay = 0, nablaz = 0, nabu, nabv, nabw ;
        Typ tmp, Fkp ;
        nablax = (phi[neib[id_f+1]]-phi[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        if((int)items[IDX_Q]==27){
        nablay = (phi[neib[id_f+2]]-phi[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;}
        nablaz = (phi[neib[id_f+up]]-phi[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;
        nabu = (phi[id_rho]*velx[id_rho]-phiold[id_rho]*uold[id_rho])/items[IDX_dt] ;
        nabv = (phi[id_rho]*vely[id_rho]-phiold[id_rho]*vold[id_rho])/items[IDX_dt] ;
        nabw = (phi[id_rho]*velz[id_rho]-phiold[id_rho]*wold[id_rho])/items[IDX_dt] ;
        phiold[id_rho] = phi[id_rho] ; uold[id_rho] = velx[id_rho] ; vold[id_rho] = vely[id_rho] ; wold[id_rho] = velz[id_rho] ;
        Typ gkx=0, gkz=0, Fk[9] ;
        for(k=0;k<items[IDX_Q];k++){
            tmp = (items[IDX_cx(k)]*velx[id_rho]  + items[IDX_cy(k)]*vely[id_rho] +  items[IDX_cz(k)]*velz[id_rho])/(items[IDX_c]*items[IDX_c]) ;

            geq [id_f+k]  = items[IDX_w(k)]*phi[id_rho] * (1.0 + 3.0*tmp + 4.5*tmp*tmp - 1.5*(velx[id_rho]*velx[id_rho]+vely[id_rho]*vely[id_rho]+velz[id_rho]*velz[id_rho])/(items[IDX_c]*items[IDX_c])) ;
            Fkp = items[IDX_dt] * items[IDX_w(k)] * (
            (1.0 - 4.0*powf(phi[id_rho]-0.5,2))/items[IDX_PFthick]*(items[IDX_cx(k)]*nablax + items[IDX_cy(k)]*nablay + items[IDX_cz(k)]*nablaz)
            / sqrtf(powf(nablax,2) + powf(nablay,2) + powf(nablaz,2) + powf(10,-8)) 
            + 3 * (items[IDX_cx(k)]*nabu+items[IDX_cy(k)]*nabv+items[IDX_cz(k)]*nabw) / powf(items[IDX_c],2) 
            ) ;
            // g[id_f+k] = g[id_f+k] - (g[id_f+k]-geq[id_f+k]+Fkp/2.0)/(0.5+0.00004*items[IDX_dt]/powf(items[IDX_dz],2)) + Fkp ; // */
            g[id_f+k] = g[id_f+k] - (g[id_f+k]-geq[id_f+k]-Fkp/1.0*0.08)/(0.5+0.08) + Fkp*0 ; // */

            gtmp[id_f+k] = g[id_f+k] ; geq[id_f+k] = 0 ; g[id_f+k] = 0 ;
        }
    }
}

template<typename Typ>
__global__ void update_scalar(Typ *items, Typ *g, Typ *scalar){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; 
    if(id_rho<items[IDX_num_calc]){
        scalar[id_rho] = 0 ; 
        for(k=0;k<items[IDX_Q];k++){
            scalar[id_rho] += g[id_f+k] ;
        }
    }
}

template<typename Typ>
__global__ void update_rho(Typ *items, Typ rhoL, Typ rhoH, Typ *f, Typ *Fx, Typ *Fy, Typ *Fz, Typ *pressure, Typ *sal, Typ *phi, Typ *rho, Typ *velx, Typ *vely, Typ *velz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    int k ; 
    if(id_rho<items[IDX_num_calc]){
        pressure[id_rho]=0 ; velx[id_rho] = 0 ; vely[id_rho] = 0 ; velz[id_rho] = 0 ; 
        for(k=0;k<items[IDX_Q];k++){
            pressure[id_rho] += f[id_f+k] ;
            velx[id_rho]+=items[IDX_cx(k)]*f[id_f+k] ;
            vely[id_rho]+=items[IDX_cy(k)]*f[id_f+k] ;
            velz[id_rho]+=items[IDX_cz(k)]*f[id_f+k] ;
        }
        rho[id_rho] = rhoL + (rhoH-rhoL)*phi[id_rho] + sal[id_rho]*0.824493 ;
        velx[id_rho] += Fx[id_rho]*items[IDX_dt]  / (2.0*rho[id_rho]) ;
        vely[id_rho] += Fy[id_rho]*items[IDX_dt]  / (2.0*rho[id_rho]) ;
        velz[id_rho] += Fz[id_rho]*items[IDX_dt]  / (2.0*rho[id_rho]) ; 
        // Fx[id_rho]=0 ; Fy[id_rho]=0 ; Fz[id_rho]=0 ; // */
    }
}
template<typename Typ>
__global__ void resetF(Typ *items, Typ *Fx, Typ *Fy, Typ *Fz, int num_F){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<num_F){
        Fx[id_rho] = 0 ; Fy[id_rho] = 0 ; Fz[id_rho] = 0 ;
    }
}

template __global__ void equ_f<float>(float*, float*, float*, float*, float*, float*) ;
template __global__ void equ_f<double>(double*, double*, double*, double*, double*, double*) ;
template __global__ void equ_g<float>(float*, float*, float*, float*, float*, float*) ;
template __global__ void equ_g<double>(double*, double*, double*, double*, double*, double*) ;
template __global__ void wall_function<float>(float*, float*, float*, int, int, int, int, int*, float*, float*, float*, float*, float*, float*) ;
template __global__ void wall_function<double>(double*, double*, double*, int, int, int, int, int*, double*, double*, double*, double*, double*, double*) ;
template __global__ void Force<float>(float*, bool, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*) ;
template __global__ void Force<double>(double*, bool, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) ;
template __global__ void col_f_MRT<float>(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*) ;
template __global__ void col_f_MRT<double>(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) ;
template __global__ void col_f_SRT<float>(float*, float*, float*, float*, float*, float*) ;
template __global__ void col_f_SRT<double>(double*, double*, double*, double*, double*, double*) ;
template __global__ void col_g_reg<float>(float*, float*, float*, float*, float*, float*, float*, float*, float*) ;
template __global__ void col_g_reg<double>(double*, double*, double*, double*, double*, double*, double*, double*, double*) ;
template __global__ void col_PF<float>(float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*) ;
template __global__ void col_PF<double>(double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) ;
template __global__ void update_scalar<float>(float*, float*, float*) ;
template __global__ void update_scalar<double>(double*, double*, double*) ;
template __global__ void update_rho<float>(float*, float, float, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*) ;
template __global__ void update_rho<double>(double*, double, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) ;
template __global__ void resetF<float>(float*, float*, float*, float*,int) ;
template __global__ void resetF<double>(double*, double*, double*, double*,int) ;
