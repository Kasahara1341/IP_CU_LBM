#include "../all.hpp"

template<typename Typ>
void IB_csv(int loop, vector<Typ>& items, vector<Typ>& posx, vector<Typ>& posy, vector<Typ>& posz, vector<Typ>& Gwx, vector<Typ>& Gwy, vector<Typ>& Gwz){
    char filename[100] ;
    sprintf(filename, "IB_point/ibm_%04d.csv", loop);
    std::ofstream file(filename);    
    file << "x,y,z,ax,ay,az\n";
    for (int i = 0 ; i<posx.size(); i++){
        file << posx[i] << "," << posy[i] << "," << posz[i] << ","
             << Gwx[i] << "," << Gwy[i] << "," << Gwz[i] << "\n";
    }
}

__device__ float profile_s(float dz, float Radius, float posx, float posy, float centerx, float centery){
    float r = Radius - 0.5*dz - sqrtf(powf(posx-centerx,2)+powf(posy-centery,2)) ;
    float result = (0 * (r<-0.5*dz)) + (0.5*(sin(3.14159*r/dz) + 1.0) * (fabsf(r)<=0.5*dz)) + (1.0 * (r>0.5*dz)) ;
    return result ;
}

template<typename Typ>
__global__ void SPM(Typ *items, Typ Radius, Typ centerx, Typ centery, Typ *f, Typ *ftmp, Typ *tau, Typ *posx, Typ *posy, Typ *Fx, Typ *Fy, Typ *Fz, Typ *velx, Typ *vely, Typ *velz){
    // smoothed-profile method
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho * (int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        Fx[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],centerx,centery) * (0 - velx[id_rho])/items[IDX_dt] ;
        Fy[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],centerx,centery) * (0 - vely[id_rho])/items[IDX_dt] ;
        Fz[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],centerx,centery) * (0 - velz[id_rho])/items[IDX_dt] ;
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] = ftmp[id_f+k] + items[IDX_w(k)]*items[IDX_dt] * 3.0
            *( items[IDX_cx(k)]*Fx[id_rho] + items[IDX_cy(k)]*Fy[id_rho] + items[IDX_cz(k)]*Fz[id_rho] )/(powf(items[IDX_c],2)) ;
            ftmp[id_f+k] = f[id_f+k] ;
        }
        // velx[id_rho] += Fx[id_rho] * items[IDX_dt] ; vely[id_rho] += Fy[id_rho] * items[IDX_dt] ; velz[id_rho] += Fz[id_rho] * items[IDX_dt] ;
        Fx[id_rho]=0 ; Fy[id_rho]=0 ; Fz[id_rho]=0 ;
    }

}

// Suzuki & Yoshino (2018) method 
template<typename Typ>
__global__ void get_IBMGw2(Typ *items, int *lattice_id, int *neib, Typ *f, Typ *tau, Typ *posx, Typ *posy, Typ *posz, Typ *posw, Typ centerx, Typ centery, Typ centerz, Typ *velx, Typ *vely, Typ *velz, Typ *Fx, Typ *Fy, Typ *Fz, Typ *Gw){
    // https://www2.nagare.or.jp/cfd/cfd32/cfd32papers/paper/F10-1.pdf 
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
        int idf[9], num_point =3, k_inv[27] ;   set_k_inv(k_inv,(int)items[IDX_Q]) ;
        Typ coeffx[5], coeffy[5], coeffz[5], posl_x[5] ;
        Typ fw[27]={0}, fw1[27]={0}, fw2[27]={0} ; 
        Typ u_w=0, v_w=0, w_w=0 ;
        idf[4]=lattice_id[id_rho] ;            
        idf[6]=neib[idf[4]*(int)items[IDX_Q]+ 9] ; idf[7]=neib[idf[4]*(int)items[IDX_Q]+2] ; idf[8]=neib[idf[4]*(int)items[IDX_Q]+7] ;
        idf[3]=neib[idf[4]*(int)items[IDX_Q]+ 3] ; idf[4] = lattice_id[id_rho]             ; idf[5]=neib[idf[4]*(int)items[IDX_Q]+1] ; 
        idf[0]=neib[idf[4]*(int)items[IDX_Q]+10] ; idf[1]=neib[idf[4]*(int)items[IDX_Q]+4] ; idf[2]=neib[idf[4]*(int)items[IDX_Q]+8] ; 
        posl_x[0] = posx[idf[3]] ; posl_x[1] = posx[idf[4]] ; posl_x[2] = posx[idf[5]] ;
        Lagrange_interpolation<Typ>(coeffx, posl_x, posw[id_rho*3+0], num_point) ;
        posl_x[0] = posy[idf[1]] ; posl_x[1] = posy[idf[4]] ; posl_x[2] = posy[idf[7]] ;
        Lagrange_interpolation<Typ>(coeffy, posl_x, posw[id_rho*3+1], num_point) ; 
        posl_x[0] = posy[idf[1]] ; posl_x[1] = posy[idf[4]] ; posl_x[2] = posy[idf[7]] ;
        Lagrange_interpolation<Typ>(coeffy, posl_x, posw[id_rho*3+1], num_point) ; 

        for(int i=0; i<num_point; i++){
            for(int j=0; j<num_point; j++){
                u_w   += coeffx[i]*coeffy[j] * velx[idf[i*num_point+j]] ;
                v_w   += coeffx[i]*coeffy[j] * vely[idf[i*num_point+j]] ;
                w_w   += coeffx[i]*coeffy[j] * velz[idf[i*num_point+j]] ;
                for(int k=0; k<items[IDX_Q];k++){
                    fw[k] += coeffx[i]*coeffy[j] * f[idf[i*num_point+j]*(int)items[IDX_Q] + k] ;
                }
            }
        }
        u_w=0 ; v_w=0 ; w_w=0 ;

        Typ unit_vectorx = (posw[id_rho*3+0]-centerx)/sqrtf( powf(posw[id_rho*3+0]-centerx,2) + powf(poswy[id_rho*3+1]-centery,2) + powf(poswz[id_rho*3+2]-centerz,2) ) ;
        Typ unit_vectory = (posw[id_rho*3+1]-centery)/sqrtf( powf(posw[id_rho*3+0]-centerx,2) + powf(poswy[id_rho*3+1]-centery,2) + powf(poswz[id_rho*3+2]-centerz,2) ) ;
        Typ unit_vectorz = (posw[id_rho*3+2]-centerz)/sqrtf( powf(posw[id_rho*3+0]-centerx,2) + powf(poswy[id_rho*3+1]-centery,2) + powf(poswz[id_rho*3+2]-centerz,2) ) ;
        Typ alpha = 22.5/180.0*3.14159, beta ;
        for(int k=0;k<items[IDX_Q];k++){
            Typ tmp  = unit_vectorx*items[IDX_cx(k)]+unit_vectory*items[IDX_cy(k)] ;
            Typ tmp2 = sqrt( powf(items[IDX_cx(k)],2) + powf(items[IDX_cy(k)],2) + powf(items[IDX_cz(k)],2) ) * sin(alpha) ;
            beta = fabs(tmp)/tmp2 ; 
            tmp2 =  0 ;
            fw1[k] = (fw[k] * (tmp<0))
            + ( ((1.0-beta)*fw[k] + beta*fw[k_inv[k]] + beta*6.0*items[IDX_w(k)]*(items[IDX_cx(k)]*u_w+items[IDX_cy(k)]*v_w+items[IDX_cz(k)]*w_w)/powf(items[IDX_c],2) ) * (0 < tmp && tmp <= tmp2) )
            + ( (fw[k_inv[k]] + 6.0*items[IDX_w(k)]*(items[IDX_cx(k)]*u_w+items[IDX_cy(k)]*v_w+items[IDX_cz(k)]*w_w)/powf(items[IDX_c],2) ) * (tmp > tmp2 )) ;
            fw2[k] = (fw[k] * (tmp>=0))
            + ( ((1.0-beta)*fw[k] + beta*fw[k_inv[k]] + beta*6.0*items[IDX_w(k)]*(items[IDX_cx(k)]*u_w+items[IDX_cy(k)]*v_w+items[IDX_cz(k)]*w_w)/powf(items[IDX_c],2)) * (-tmp2 < tmp && tmp < 0))
            + ( (fw[k_inv[k]] + 6.0*items[IDX_w(k)]*(items[IDX_cx(k)]*u_w+items[IDX_cy(k)]*v_w+items[IDX_cz(k)]*w_w)/powf(items[IDX_c],2) ) * (tmp < - tmp2)) ;
        }
        fw1[0]=fw[0] ; fw2[0]=fw[0] ; fw1[5]=fw[5] ; fw2[5]=fw[5] ; fw1[6]=fw[6] ; fw2[6]=fw[6] ; 

        Typ sigma_x=0, sigma_y=0, sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma x
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*u_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *1  - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw2[k]-items[IDX_w(k)])*u_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*v_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0  - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw2[k]-items[IDX_w(k)])*u_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*w_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0  - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw2[k]-items[IDX_w(k)])*u_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+0] = -(sigma_x*unit_vectorx + sigma_y*unit_vectory)/items[IDX_dz] ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma y
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*u_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw2[k]-items[IDX_w(k)])*v_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*v_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *1 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw2[k]-items[IDX_w(k)])*v_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*w_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw2[k]-items[IDX_w(k)])*v_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+1] = -(sigma_x*unit_vectorx + sigma_y*unit_vectory)/items[IDX_dz] ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma z
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*u_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw2[k]-items[IDX_w(k)])*w_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*v_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *0 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw2[k]-items[IDX_w(k)])*w_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*w_w)/powf(items[IDX_c],2)
                     - 0.5/tau[lattice_id[id_rho]] * fw2[k]/3.0 *1 - (1.0-0.5/tau[lattice_id[id_rho]]) * (fw2[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw2[k]-items[IDX_w(k)])*w_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gwz[id_rho*3+2] = -(sigma_x*unit_vectorx + sigma_y*unit_vectory)/items[IDX_dz] ;

        for(int i=0; i<num_point; i++){
            for(int j=0; j<num_point; j++){
                atomicAdd(&Fx[idf[i*num_point+j]], Gwx[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],2)*coeffx[i]*coeffy[j]) ;
                atomicAdd(&Fy[idf[i*num_point+j]], Gwy[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],2)*coeffx[i]*coeffy[j]) ;
                atomicAdd(&Fz[idf[i*num_point+j]], Gwz[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],2)*coeffx[i]*coeffy[j]) ;
            }
        }
        sigma_x=0, sigma_y=0, sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma x
            sigma_x += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*u_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+0] = sigma_x*unit_vectorx + sigma_y*unit_vectory + sigma_z*unit_vectorz ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma y
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*v_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+1] = sigma_x*unit_vectorx + sigma_y*unit_vectory + sigma_z*unit_vectorz ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma z
            sigma_x += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cx(k)]-u_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cy(k)]-v_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z += 0.5/tau[lattice_id[id_rho]] * fw1[k]/3.0 *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * (fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cz(k)]-w_w)-(1.0*fw1[k]-items[IDX_w(k)])*w_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+2] = sigma_x*unit_vectorx + sigma_y*unit_vectory + sigma_z*unit_vectorz ;
    }
}

__global__ void moving_body(float *items, int *IB_index, float *FB, float *Torque, float *velB, float *Gwx, float *Gwy, float *Gwz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
    }
}

// direct forcing method
template<typename Typ>
__global__ void get_IBMGw(Typ *items, int *lattice_id, int *neib, Typ *velx, Typ *vely, Typ *velz, Typ *posx, Typ *posy, Typ *posz, Typ *poswx, Typ *poswy, Typ *poswz, Typ *Fx, Typ *Fy, Typ *Fz, Typ *Gwx, Typ *Gwy, Typ *Gwz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
        int idf[25], num_point =3 ;  
        Typ coeffx[5], coeffy[5], posl_x[5] ;
        Typ u_w=0, v_w=0, w_w=0 ; 
        if(num_point==5){
            idf[12]=lattice_id[id_rho] ;
            idf[13] = neib[idf[12]*(int)items[IDX_Q] +  1] ; idf[14] = neib[idf[13]*(int)items[IDX_Q] +  1] ;
            idf[11] = neib[idf[12]*(int)items[IDX_Q] +  3] ; idf[10] = neib[idf[11]*(int)items[IDX_Q] +  3] ;
            idf[17] = neib[idf[12]*(int)items[IDX_Q] +  2] ; idf[22] = neib[idf[17]*(int)items[IDX_Q] +  2] ;
            idf[ 7] = neib[idf[12]*(int)items[IDX_Q] +  4] ; idf[ 2] = neib[idf[ 7]*(int)items[IDX_Q] +  4] ;
            idf[18] = neib[idf[12]*(int)items[IDX_Q] +  7] ; idf[24] = neib[idf[18]*(int)items[IDX_Q] +  7] ;
            idf[ 6] = neib[idf[12]*(int)items[IDX_Q] + 10] ; idf[ 0] = neib[idf[ 6]*(int)items[IDX_Q] + 10] ;
            idf[ 8] = neib[idf[12]*(int)items[IDX_Q] +  8] ; idf[ 4] = neib[idf[ 8]*(int)items[IDX_Q] +  8] ;
            idf[16] = neib[idf[12]*(int)items[IDX_Q] +  9] ; idf[20] = neib[idf[16]*(int)items[IDX_Q] +  9] ;
            idf[19] = neib[idf[14]*(int)items[IDX_Q] +  2] ; idf[ 9] = neib[idf[14]*(int)items[IDX_Q] +  4] ;
            idf[15] = neib[idf[10]*(int)items[IDX_Q] +  2] ; idf[ 5] = neib[idf[10]*(int)items[IDX_Q] +  4] ;
            idf[23] = neib[idf[22]*(int)items[IDX_Q] +  1] ; idf[21] = neib[idf[22]*(int)items[IDX_Q] +  3] ;
            idf[ 3] = neib[idf[ 2]*(int)items[IDX_Q] +  1] ; idf[ 1] = neib[idf[ 2]*(int)items[IDX_Q] +  3] ;    
            posl_x[0] = posx[idf[10]] ; posl_x[1] = posx[idf[11]] ; posl_x[2] = posx[idf[12]] ;
            posl_x[3] = posx[idf[13]] ; posl_x[4] = posx[idf[14]] ;
            Lagrange_interpolation<Typ>(coeffx, posl_x, poswx[id_rho], 5) ;
            posl_x[0] = posy[idf[ 2]] ; posl_x[1] = posy[idf[ 7]] ; posl_x[2] = posy[idf[12]] ;
            posl_x[3] = posy[idf[17]] ; posl_x[4] = posy[idf[22]] ;
            Lagrange_interpolation<Typ>(coeffy, posl_x, poswy[id_rho], 5) ; 
        } 
        else if(num_point==3){
            idf[4]=lattice_id[id_rho] ;            
            idf[6]=neib[idf[4]*(int)items[IDX_Q]+ 9] ; idf[7]=neib[idf[4]*(int)items[IDX_Q]+2] ; idf[8]=neib[idf[4]*(int)items[IDX_Q]+7] ;
            idf[3]=neib[idf[4]*(int)items[IDX_Q]+ 3] ; idf[4] = lattice_id[id_rho]             ; idf[5]=neib[idf[4]*(int)items[IDX_Q]+1] ; 
            idf[0]=neib[idf[4]*(int)items[IDX_Q]+10] ; idf[1]=neib[idf[4]*(int)items[IDX_Q]+4] ; idf[2]=neib[idf[4]*(int)items[IDX_Q]+8] ; 
            posl_x[0] = posx[idf[3]] ; posl_x[1] = posx[idf[4]] ; posl_x[2] = posx[idf[5]] ;
            Lagrange_interpolation<Typ>(coeffx, posl_x, poswx[id_rho], num_point) ;
            posl_x[0] = posy[idf[1]] ; posl_x[1] = posy[idf[4]] ; posl_x[2] = posy[idf[7]] ;
            Lagrange_interpolation<Typ>(coeffy, posl_x, poswy[id_rho], num_point) ; 
        }
        for(int i=0; i<num_point; i++){
            for(int j=0; j<num_point; j++){
                u_w += coeffx[i]*coeffy[j] * velx[idf[i*num_point+j]] ;
                v_w += coeffx[i]*coeffy[j] * vely[idf[i*num_point+j]] ;
                w_w += coeffx[i]*coeffy[j] * velz[idf[i*num_point+j]] ;
            }
        }
        Gwx[id_rho] += - u_w / items[IDX_dt] ; Gwy[id_rho] += - v_w / items[IDX_dt] ; Gwz[id_rho] += - w_w / items[IDX_dt] ;
        for(int i=0; i<num_point; i++){
            for(int j=0; j<num_point; j++){
                atomicAdd(&Fx[idf[i*num_point+j]], Gwx[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],1)*coeffx[i]*coeffy[j]) ;
                atomicAdd(&Fy[idf[i*num_point+j]], Gwy[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],1)*coeffx[i]*coeffy[j]) ;
                atomicAdd(&Fz[idf[i*num_point+j]], Gwz[id_rho]*powf(items[IDX_dIBM]/items[IDX_dz],1)*coeffx[i]*coeffy[j]) ;
            }
        }
    }
}

template<typename Typ>
__global__ void update_velIBM(Typ *items, int *lattice_id, Typ *f, Typ *ftmp, Typ *pressure, Typ *tau, Typ *velx, Typ *vely, Typ *velz, Typ *velx_old, Typ *vely_old, Typ *velz_old, Typ *Fx, Typ *Fy, Typ *Fz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_calc]){
        int id_f = id_rho * (int)items[IDX_Q] ;
        pressure[id_rho]=0 ; velx_old[id_rho] = 0 ; vely_old[id_rho] = 0 ; velz_old[id_rho] = 0 ; 
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] = ftmp[id_f+k] + (1.0-0.5/tau[id_rho])* items[IDX_w(k)]*items[IDX_dt]*
            3.0*(items[IDX_cx(k)]*Fx[id_rho]+items[IDX_cy(k)]*Fy[id_rho]+items[IDX_cy(k)]*Fy[id_rho])/(powf(items[IDX_c],2)) ;
            pressure[id_rho] += f[id_f+k] ;
            velx_old[id_rho]+=items[IDX_cx(k)]*f[id_f+k] ;
            vely_old[id_rho]+=items[IDX_cy(k)]*f[id_f+k] ;
            velz_old[id_rho]+=items[IDX_cz(k)]*f[id_f+k] ;
        } // */
        velx[id_rho] = velx_old[id_rho] + items[IDX_dt] * Fx[id_rho]/2.0  ;
        vely[id_rho] = vely_old[id_rho] + items[IDX_dt] * Fy[id_rho]/2.0  ;
        velz[id_rho] = velz_old[id_rho] + items[IDX_dt] * Fz[id_rho]/2.0  ;
        // Fx[id_rho]=0 ; Fy[id_rho]=0 ; Fz[id_rho]=0;
    }
}


template void IB_csv<float>(int, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&);
template void IB_csv<double>(int,vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&);
template __global__ void SPM<float>(float*, float, float, float, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void SPM<double>(double*, double, double, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
template __global__ void get_IBMGw2<float>(float*, int*, int*, float*, float*, float*, float*, float*, float*, float, float, float, float*, float*, float*, float*, float*, float*, float*);
template __global__ void get_IBMGw2<double>(double*, int*, int*, double*, double*, double*, double*, double*, double*, double, double, double, double*, double*, double*, double*, double*, double*, double*);
template __global__ void get_IBMGw<float>(float*, int*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void get_IBMGw<double>(double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
template __global__ void update_velIBM<float>(float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void update_velIBM<double>(double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
