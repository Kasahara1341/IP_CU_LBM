#include "../all.hpp"

template<typename Typ>
void IB_csv(int loop, vector<Typ>& items, vector<Typ>& pos, vector<Typ>& velw, vector<Typ>& Gw){
    char filename[100] ;
    sprintf(filename, "IB_point/ibm_%04d.csv", loop);
    std::ofstream file(filename);    
    file << "x,y,z,velx,vely,velz,ax,ay,az\n";
    for (int i = 0 ; i<pos.size()/3; i++){
        file << pos[i*3] << "," << pos[i*3+1] << "," << pos[i*3+2] << ","
             << velw[i*3] << "," << velw[i*3+1] << "," << velw[i*3+2] << ","
             << Gw[i*3] << "," << Gw[i*3+1] << "," << Gw[i*3+2] << "\n";
    }
}

__host__ __device__ void set_quaternionS(int IB_index, float q0, float q1, float q2, float q3, float *quaS){
            quaS[IB_index*9+0] = powf(q0,2)-powf(q2,2)-powf(q3,2)+powf(q1,2) ;
            quaS[IB_index*9+1] = 2*(q1*q2 - q0*q3) ; 
            quaS[IB_index*9+2] = 2*(q1*q3 + q0*q2) ;

            quaS[IB_index*9+3] = 2*(q1*q2 + q0*q3) ; 
            quaS[IB_index*9+4] = powf(q0,2)-powf(q3,2)-powf(q1,2)+powf(q2,2) ;
            quaS[IB_index*9+5] = 2*(q2*q3 - q0*q1) ;

            quaS[IB_index*9+6] = 2*(q1*q3 - q0*q2) ; 
            quaS[IB_index*9+7] = 2*(q2*q3 + q0*q1) ;
            quaS[IB_index*9+8] = powf(q0,2)-powf(q1,2)-powf(q2,2)+powf(q3,2) ;
}
template<typename Typ>
void set_quaternionS(int IB_index, Typ q0, Typ q1, Typ q2, Typ q3, vector<Typ>& quaS){
            quaS[IB_index*9+0] = powf(q0,2)-powf(q2,2)-powf(q3,2)+powf(q1,2) ;
            quaS[IB_index*9+1] = 2*(q1*q2 - q0*q3) ; 
            quaS[IB_index*9+2] = 2*(q1*q3 + q0*q2) ;

            quaS[IB_index*9+3] = 2*(q1*q2 + q0*q3) ; 
            quaS[IB_index*9+4] = powf(q0,2)-powf(q3,2)-powf(q1,2)+powf(q2,2) ;
            quaS[IB_index*9+5] = 2*(q2*q3 - q0*q1) ;

            quaS[IB_index*9+6] = 2*(q1*q3 - q0*q2) ; 
            quaS[IB_index*9+7] = 2*(q2*q3 + q0*q1) ;
            quaS[IB_index*9+8] = powf(q0,2)-powf(q1,2)-powf(q2,2)+powf(q3,2) ;
}

__device__ float profile_s(float dz, float Radius, float posx, float posy, float centerx, float centery){
    float r = Radius - 0.5*dz - sqrtf(powf(posx-centerx,2)+powf(posy-centery,2)) ;
    float result = (0 * (r<-0.5*dz)) + (0.5*(sin(3.14159*r/dz) + 1.0) * (fabsf(r)<=0.5*dz)) + (1.0 * (r>0.5*dz)) ;
    return result ;
}

__device__ float profile_s2(float limit_lenght, float distance){
    float r = distance, in_line, out_line ; 
    in_line  = 1.0-limit_lenght + powf(limit_lenght/2,2) ;
    out_line = 1.0+limit_lenght + powf(limit_lenght/2,2) ;
    float result =    (0 * ( out_line < r)) 
                    + (0.5*(sin(3.14159*r/limit_lenght) + 1.0) * ( in_line <= r && r<= out_line)) 
                    + (1.0 * (r < in_line)) ;
    return result ;
}

__device__ float weightFunction(float dx, float posx, float posl_x){
    float r = fabsf(posx - posl_x)/dx ;
    float result = (3.0-2.0*r+sqrtf((1.0+4.0*r-4.0*r*r)*(r<1.0)))*(r<1.0)/8.0 
        + (5.0-2.0*r-sqrtf((-7.0+12.0*r-4.0*r*r)*(1.0<=r&&r<=2.0)))*(1.0<=r&&r<=2.0)/8.0 ;
    return result ;
}

template<typename Typ>
__global__ void SPM(Typ *items, Typ Radius, Typ *posB, Typ *f, Typ *ftmp, Typ *tau, Typ *posx, Typ *posy, Typ *Fx, Typ *Fy, Typ *Fz, Typ *velx, Typ *vely, Typ *velz, Typ *velw){
    // smoothed-profile method
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho * (int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        Fx[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],posB[0*3+0],posB[0*3+1]) * (velw[0] - velx[id_rho])/items[IDX_dz] ;
        Fy[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],posB[0*3+0],posB[0*3+1]) * (velw[1] - vely[id_rho])/items[IDX_dz] ;
        Fz[id_rho] = profile_s(items[IDX_dz],Radius,posx[id_rho],posy[id_rho],posB[0*3+0],posB[0*3+1]) * (velw[2] - velz[id_rho])/items[IDX_dz] ;
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] = ftmp[id_f+k] + items[IDX_w(k)]*items[IDX_dt] * 3.0
            *( items[IDX_cx(k)]*Fx[id_rho] + items[IDX_cy(k)]*Fy[id_rho] + items[IDX_cz(k)]*Fz[id_rho] )/(powf(items[IDX_c],2)) ;
            ftmp[id_f+k] = f[id_f+k] ;
        }
        // velx[id_rho] += Fx[id_rho] * items[IDX_dt] ; vely[id_rho] += Fy[id_rho] * items[IDX_dt] ; velz[id_rho] += Fz[id_rho] * items[IDX_dt] ;
        Fx[id_rho]=0 ; Fy[id_rho]=0 ; Fz[id_rho]=0 ;
    }
}

template<typename Typ>
__global__ void SPM_ellipse(Typ *items, Typ Rada, Typ Radb, Typ *quaS, Typ *posB, Typ *f, Typ *tau, Typ *posx, Typ *posy, Typ *posz, Typ *velx, Typ *vely, Typ *velz, Typ *velB, Typ *angleVB){
    // smoothed-profile method
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho * (int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        float X1, Y1, Z1 ;
        float velBx, velBy, velBz ;
        X1 = quaS[0]*(posx[id_rho]-posB[0]) + quaS[3]*(posy[id_rho]-posB[1]) + quaS[6]*(posz[id_rho]-posB[2]) ;
        Y1 = quaS[1]*(posx[id_rho]-posB[0]) + quaS[4]*(posy[id_rho]-posB[1]) + quaS[7]*(posz[id_rho]-posB[2]) ;
        Z1 = quaS[2]*(posx[id_rho]-posB[0]) + quaS[5]*(posy[id_rho]-posB[1]) + quaS[8]*(posz[id_rho]-posB[2]) ;
        velBx = velB[0] + quaS[0]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[3]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[6]*(angleVB[0]*Y1-angleVB[1]*X1) ;
        velBy = velB[1] + quaS[1]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[4]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[7]*(angleVB[0]*Y1-angleVB[1]*X1) ;
        velBz = velB[2] + quaS[2]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[5]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[8]*(angleVB[0]*Y1-angleVB[1]*X1) ;

        X1 = quaS[0]*(posx[id_rho]-posB[0]) + quaS[1]*(posy[id_rho]-posB[1]) + quaS[2]*(posz[id_rho]-posB[2]) ;
        Y1 = quaS[3]*(posx[id_rho]-posB[0]) + quaS[4]*(posy[id_rho]-posB[1]) + quaS[5]*(posz[id_rho]-posB[2]) ;
        Z1 = quaS[6]*(posx[id_rho]-posB[0]) + quaS[7]*(posy[id_rho]-posB[1]) + quaS[8]*(posz[id_rho]-posB[2]) ;
        // velBx = velB[0] + quaS[0]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[1]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[2]*(angleVB[0]*Y1-angleVB[1]*X1) ;
        // velBy = velB[1] + quaS[3]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[4]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[5]*(angleVB[0]*Y1-angleVB[1]*X1) ;
        // velBz = velB[2] + quaS[6]*(angleVB[1]*Z1-angleVB[2]*Y1) + quaS[7]*(angleVB[2]*X1-angleVB[0]*Z1) + quaS[8]*(angleVB[0]*Y1-angleVB[1]*X1) ;
        if(items[IDX_Q]==27){Z1 = Y1 ;}// 3D(xy face)
        float distance = powf(X1/Rada,2) + powf(Z1/Radb,2), fx, fy, fz ; 
        fx = profile_s2(1.0,distance) * (velBx - velx[id_rho])/items[IDX_dt] ;
        fy = profile_s2(1.0,distance) * (velBy - vely[id_rho])/items[IDX_dt] ;
        fz = profile_s2(1.0,distance) * (velBz - velz[id_rho])/items[IDX_dt] ;

        fx = profile_s2(items[IDX_dz]/Rada,distance) * (velB[0] - velx[id_rho])/items[IDX_dt] ;
        fy = profile_s2(items[IDX_dz]/Rada,distance) * (velB[1] - vely[id_rho])/items[IDX_dt] ;
        fz = profile_s2(items[IDX_dz]/Rada,distance) * (velB[2] - velz[id_rho])/items[IDX_dt] ;
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] += items[IDX_w(k)]*items[IDX_dt] * 3.0
            *( items[IDX_cx(k)]*fx + items[IDX_cy(k)]*fy + items[IDX_cz(k)]*fz )/(powf(items[IDX_c],2)) ;
        }
    }
}

// Suzuki & Yoshino (2018) method 
template<typename Typ>
__global__ void get_IBMGw2(Typ *items, int *lattice_id, int *neib, Typ *f, Typ *tau, Typ *posx, Typ *posy, Typ *posz, Typ *posw, Typ *posB, Typ *nBvec, 
    Typ *velx, Typ *vely, Typ *velz, Typ *velw, Typ *Fx, Typ *Fy, Typ *Fz, Typ *Gw, Typ rhof){
    // https://www2.nagare.or.jp/cfd/cfd32/cfd32papers/paper/F10-1.pdf 
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
        int num_point =3, k_inv[27] ;   set_k_inv(k_inv,(int)items[IDX_Q]) ;
        int idfIndex[3][3][3]={0} ; 
        if(items[IDX_Q]==9){
            idfIndex[0][1][2]=6 ; idfIndex[1][1][2]=2 ; idfIndex[2][1][2]=5 ;
            idfIndex[0][1][1]=3 ; idfIndex[1][1][1]=0 ; idfIndex[2][1][1]=1 ;
            idfIndex[0][1][0]=5 ; idfIndex[1][1][0]=4 ; idfIndex[2][1][0]=8 ;
        }
        if(items[IDX_Q]==27){
            idfIndex[0][0][0]=26 ; idfIndex[0][0][1]=10 ; idfIndex[0][0][2]=25 ; idfIndex[0][1][0]=14 ; idfIndex[0][1][1]=3 ; idfIndex[0][1][2]=13 ; idfIndex[0][2][0]=24 ; idfIndex[0][2][1]=9  ; idfIndex[0][2][2]=23 ; 
            idfIndex[1][0][0]=18 ; idfIndex[1][0][1]=4  ; idfIndex[1][0][2]=17 ; idfIndex[1][1][0]=6  ; idfIndex[1][1][1]=0 ; idfIndex[1][1][2]=5  ; idfIndex[1][2][0]=16 ; idfIndex[1][2][1]=2  ; idfIndex[1][2][2]=15 ; 
            idfIndex[2][0][0]=22 ; idfIndex[2][0][1]=8  ; idfIndex[2][0][2]=21 ; idfIndex[2][1][0]=12 ; idfIndex[2][1][1]=1 ; idfIndex[2][1][2]=11 ; idfIndex[2][2][0]=20 ; idfIndex[2][2][1]=7  ; idfIndex[2][2][2]=19 ; 
        }
        Typ coeffx[3], coeffy[3], coeffz[3] ;
        Typ fw[27]={0}, fw1[27]={0}, fw2[27]={0} ;
        Typ u_w=velw[id_rho*3+0], v_w=velw[id_rho*3+1], w_w=velw[id_rho*3+2] ;
        // coeffx[0] = weightFunction(items[IDX_dz],posw[id_rho*3+0],posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+3]]) ;
        // coeffx[1] = weightFunction(items[IDX_dz],posw[id_rho*3+0],posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]]) ;
        // coeffx[2] = weightFunction(items[IDX_dz],posw[id_rho*3+0],posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+1]]) ;
        Typ posl_x[3] ;
        posl_x[0] = posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+3]] ;
        posl_x[1] = posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]] ;
        posl_x[2] = posx[neib[lattice_id[id_rho]*(int)items[IDX_Q]+1]] ;
        Lagrange_interpolation<Typ>(coeffx, posl_x, posw[id_rho*3+0], num_point) ;
        if(items[IDX_Q]==27){
            // coeffy[0] = weightFunction(items[IDX_dz],posw[id_rho*3+1],posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+4]]) ;
            // coeffy[1] = weightFunction(items[IDX_dz],posw[id_rho*3+1],posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]]) ;
            // coeffy[2] = weightFunction(items[IDX_dz],posw[id_rho*3+1],posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+2]]) ;
            // coeffz[0] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+6]]) ;
            // coeffz[1] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]]) ;
            // coeffz[2] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+5]]) ;
            // Typ sum_coef=0;
            // for(int i=0;i<3;i++){sum_coef += coeffx[i] ;} for(int i=0;i<3;i++){coeffx[i] /= sum_coef ;}
            // sum_coef=0;
            // for(int i=0;i<3;i++){sum_coef += coeffy[i] ;} for(int i=0;i<3;i++){coeffy[i] /= sum_coef ;}
            // sum_coef=0;
            // for(int i=0;i<3;i++){sum_coef += coeffz[i] ;} for(int i=0;i<3;i++){coeffz[i] /= sum_coef ;}
            posl_x[0] = posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+4]] ; 
            posl_x[1] = posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]] ; 
            posl_x[2] = posy[neib[lattice_id[id_rho]*(int)items[IDX_Q]+2]] ;
            Lagrange_interpolation<Typ>(coeffy, posl_x, posw[id_rho*3+1], num_point) ;
            posl_x[0] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+6]] ; 
            posl_x[1] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]] ; 
            posl_x[2] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+5]] ;
            Lagrange_interpolation<Typ>(coeffz, posl_x, posw[id_rho*3+2], num_point) ; 
        }
        else{
            coeffy[0]=0;coeffy[1]=1;coeffy[2]=0;
            // coeffz[0] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+4]]) ;
            // coeffz[1] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]]) ;
            // coeffz[2] = weightFunction(items[IDX_dz],posw[id_rho*3+2],posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+2]]) ;            
            // Typ sum_coef=0;
            // for(int i=0;i<3;i++){sum_coef += coeffx[i] ;} for(int i=0;i<3;i++){coeffx[i] /= sum_coef ;}
            // sum_coef=0;
            // for(int i=0;i<3;i++){sum_coef += coeffz[i] ;} for(int i=0;i<3;i++){coeffz[i] /= sum_coef ;}
            posl_x[0] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+4]] ;
            posl_x[1] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+0]] ;
            posl_x[2] = posz[neib[lattice_id[id_rho]*(int)items[IDX_Q]+2]] ;
            Lagrange_interpolation<Typ>(coeffz, posl_x, posw[id_rho*3+2], num_point) ;
        }

        for(int i=0; i<num_point; i++){
            for(int j=0 + (items[IDX_Q]==9); j<num_point - (items[IDX_Q]==9); j++){
                for(int l=0;l<num_point;l++){
                    int id=neib[lattice_id[id_rho]*(int)items[IDX_Q]+idfIndex[i][j][l]] ;               
                    for(int k=0; k<items[IDX_Q];k++){
                        fw[k] += coeffx[i]*coeffy[j]*coeffz[l] * f[id*(int)items[IDX_Q] + k] * (id<items[IDX_num_calc]) ;
                    }
                }
            }
        }
        // int check_id = 5 ; if(id_rho==check_id){printf("k\n"); }
        Typ alpha = 22.5/180.0*3.14159, beta ;
        for(int k=0;k<items[IDX_Q];k++){
            Typ tmp  = (nBvec[id_rho*3+0]*items[IDX_cx(k)] + nBvec[id_rho*3+1]*items[IDX_cy(k)] + nBvec[id_rho*3+2]*items[IDX_cz(k)])/items[IDX_c] ; // c times n
            Typ tmp2 = sqrtf( powf(items[IDX_cx(k)],2) + powf(items[IDX_cy(k)],2) + powf(items[IDX_cz(k)],2))/items[IDX_c] * sin(alpha) ; // |c|sin(alpha)
            Typ tmp3 = 6.0*items[IDX_w(k)]* (items[IDX_cx(k)]*u_w + items[IDX_cy(k)]*v_w + items[IDX_cz(k)]*w_w)/powf(items[IDX_c],2) ;   // 6EcU
            beta     = fabs(tmp)/(tmp2+powf(10,-10)) ; 
            fw1[k] = (fw[k] * (tmp<=0))
            + ( ((1.0-beta)*fw[k] + beta*(fw[k_inv[k]]+tmp3) ) * (0 < tmp && tmp <= tmp2) )
            + ( (fw[k_inv[k]] + tmp3 ) * (tmp > tmp2 )) ;
            fw2[k] = (fw[k] * (tmp>=0))
            + ( ((1.0-beta)*fw[k] + beta*(fw[k_inv[k]]+tmp3) ) * (-tmp2 <= tmp && tmp < 0))
            + ( (fw[k_inv[k]] + tmp3 ) * (tmp < - tmp2)) ;
            // if(id_rho==check_id){
                // printf("k=%d fw=%e fw1=%e fw2=%e  tmp=%e beta=%e |c|sin(alpha)=%f\n",k,fw[k],fw1[k],fw2[k],tmp,beta,tmp2); 
                // printf("fw1 [%d]  fw2 [%d] \n",(0 < tmp && tmp <= tmp2),(-tmp2 <= tmp && tmp < 0)) ;
            // }
        }


        Typ sigma_x=0, sigma_y=0, sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma x
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *1  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cx(k)]-u_w)*(items[IDX_cx(k)]-u_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*u_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cx(k)]-u_w)*(items[IDX_cy(k)]-v_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*u_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cx(k)]-u_w)*(items[IDX_cz(k)]-w_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*u_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+0] = -(sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2]) ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma y
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cy(k)]-v_w)*(items[IDX_cx(k)]-u_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*v_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *1  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cy(k)]-v_w)*(items[IDX_cy(k)]-v_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*v_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cy(k)]-v_w)*(items[IDX_cz(k)]-w_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*v_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+1] = -(sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2]) ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma z
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cz(k)]-w_w)*(items[IDX_cx(k)]-u_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*w_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *0  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cz(k)]-w_w)*(items[IDX_cy(k)]-v_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*w_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * (fw1[k]-fw2[k]) *1  + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*((fw1[k]-fw2[k])*(items[IDX_cz(k)]-w_w)*(items[IDX_cz(k)]-w_w)-(3*(fw1[k]-fw2[k])-items[IDX_w(k)])*w_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+2] = -(sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2])  ;


        int dimension=3 ; if(items[IDX_Q]==9){dimension=2;}
        for(int i=0; i<num_point; i++){
            for(int j=0 + (items[IDX_Q]==9); j < num_point - (items[IDX_Q]==9) ; j++){
                for(int l=0;l<num_point;l++){
                    int id=neib[lattice_id[id_rho]*(int)items[IDX_Q]+idfIndex[i][j][l]] ;
                    atomicAdd(&Fx[id], Gw[id_rho*3+0]*powf(items[IDX_dIBM]/items[IDX_dz],dimension-1)/items[IDX_dz]*coeffx[i]*coeffy[j]*coeffz[l] *(id<items[IDX_num_calc]) ) ;
                    atomicAdd(&Fy[id], Gw[id_rho*3+1]*powf(items[IDX_dIBM]/items[IDX_dz],dimension-1)/items[IDX_dz]*coeffx[i]*coeffy[j]*coeffz[l] *(id<items[IDX_num_calc]) ) ;
                    atomicAdd(&Fz[id], Gw[id_rho*3+2]*powf(items[IDX_dIBM]/items[IDX_dz],dimension-1)/items[IDX_dz]*coeffx[i]*coeffy[j]*coeffz[l] *(id<items[IDX_num_calc]) ) ;
                }
            }
        }
        sigma_x=0, sigma_y=0, sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma x
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cx(k)]-u_w)-(3*fw1[k]-items[IDX_w(k)])*u_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cy(k)]-v_w)-(3*fw1[k]-items[IDX_w(k)])*u_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cx(k)]-u_w)*(items[IDX_cz(k)]-w_w)-(3*fw1[k]-items[IDX_w(k)])*u_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+0] = sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2] ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma y
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cx(k)]-u_w)-(3*fw1[k]-items[IDX_w(k)])*v_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cy(k)]-v_w)-(3*fw1[k]-items[IDX_w(k)])*v_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cy(k)]-v_w)*(items[IDX_cz(k)]-w_w)-(3*fw1[k]-items[IDX_w(k)])*v_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+1] = sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2] ;
        sigma_x=0 ; sigma_y=0 ; sigma_z=0 ;
        for(int k=0;k<items[IDX_Q];k++){ // sigma z
            sigma_x -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cx(k)]-u_w)-(3*fw1[k]-items[IDX_w(k)])*w_w*u_w)/powf(items[IDX_c],2) ;
            sigma_y -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *0 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cy(k)]-v_w)-(3*fw1[k]-items[IDX_w(k)])*w_w*v_w)/powf(items[IDX_c],2) ;
            sigma_z -= 0.5/tau[lattice_id[id_rho]] * fw1[k] *1 + (1.0-0.5/tau[lattice_id[id_rho]]) * 3*(fw1[k]*(items[IDX_cz(k)]-w_w)*(items[IDX_cz(k)]-w_w)-(3*fw1[k]-items[IDX_w(k)])*w_w*w_w)/powf(items[IDX_c],2) ;
        }
        Gw[id_rho*3+2] = sigma_x*nBvec[id_rho*3+0] + sigma_y*nBvec[id_rho*3+1] + sigma_z*nBvec[id_rho*3+2] ;
        for(int i=0;i<3;i++){
            Gw[id_rho*3+i]*= rhof*powf(items[IDX_c],2)/3.0 * powf(items[IDX_dIBM],dimension-1) ;
        }
    }
}

template<typename Typ>
__global__ void update_velIBM(Typ *items, int *lattice_id, Typ *f, Typ *ftmp, Typ *pressure, Typ *tau, Typ *velx, Typ *vely, Typ *velz, Typ *velx_old, Typ *vely_old, Typ *velz_old, Typ *Fx, Typ *Fy, Typ *Fz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_calc]){
        int id_f = id_rho * (int)items[IDX_Q] ;
        pressure[id_rho]=0 ; velx[id_rho] = 0 ; vely[id_rho] = 0 ; velz[id_rho] = 0 ; 
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] = f[id_f+k] + (1.0-0.5/tau[id_rho]*0)* items[IDX_w(k)]*items[IDX_dt]*
                (items[IDX_cx(k)]*Fx[id_rho] + items[IDX_cy(k)]*Fy[id_rho] + items[IDX_cz(k)]*Fz[id_rho])  ;
            pressure[id_rho]+= f[id_f+k] ;
            velx[id_rho]+=items[IDX_cx(k)]*f[id_f+k] ;
            vely[id_rho]+=items[IDX_cy(k)]*f[id_f+k] ;
            velz[id_rho]+=items[IDX_cz(k)]*f[id_f+k] ;
        } // */
        // velx[id_rho] += powf(items[IDX_c],2)*items[IDX_dt] * Fx[id_rho]/2.0 ;
        // vely[id_rho] += powf(items[IDX_c],2)*items[IDX_dt] * Fy[id_rho]/2.0 ;
        // velz[id_rho] += powf(items[IDX_c],2)*items[IDX_dt] * Fz[id_rho]/2.0 ;
        // Fx[id_rho]=0 ; Fy[id_rho]=0 ; Fz[id_rho]=0;
    }
}

__global__ void update_IBbody(float *items, int IB_index, float *massB, float *densB, float *inertia, float *FB, float *posB, float *Torque, float *velB, float *quat, 
    float *quaS, float *angleVB, float *posw, float *Gw, float *quatold, float rhof){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<1){
        float quaSold[9], FBold[3], Torqueold[3], velBold[3], angleVBold[3] ;
        set_quaternionS(0,quatold[IB_index*4+0],quatold[IB_index*4+1],
            quatold[IB_index*4+2],quatold[IB_index*4+3],quaSold) ;
        for(int i=0;i<3;i++){
            FBold[i] = FB[IB_index*3+i] ; FB[IB_index*3+i]=0 ;
            Torqueold[i] = Torque[IB_index*3+i] ; Torque[IB_index*3+i]=0 ;
            velBold[i] = velB[IB_index*3+i] ; 
            angleVBold[i] = angleVB[IB_index*3+i] ;
        }

        // update Force & Torque
        for(int i=0;i<items[IDX_num_IBPoints];i++){
            FB[IB_index*3+0] += Gw[i*3+0] ;
            FB[IB_index*3+1] += Gw[i*3+1] ;
            FB[IB_index*3+2] += Gw[i*3+2] ;
            Torque[IB_index*3+0] += (posw[i*3+1]-posB[IB_index*3+1])*Gw[i*3+2] - (posw[i*3+2]-posB[IB_index*3+2])*Gw[i*3+1] ;
            Torque[IB_index*3+1] += (posw[i*3+2]-posB[IB_index*3+2])*Gw[i*3+0] - (posw[i*3+0]-posB[IB_index*3+0])*Gw[i*3+2] ;
            Torque[IB_index*3+2] += (posw[i*3+0]-posB[IB_index*3+0])*Gw[i*3+1] - (posw[i*3+1]-posB[IB_index*3+1])*Gw[i*3+0] ;
        }
        FB[IB_index*3+0] += (1.0-rhof/densB[IB_index])*massB[IB_index]*9.81 ;
        // FB[IB_index*3+2] += (1.0-rhof/densB[IB_index])*massB[IB_index]*9.81 ;
        
        // Finalyze FB & update velocity of IB_body
        for(int i=0;i<3;i++){
            velB[IB_index*3+i] += items[IDX_dt]*(3.0*FB[IB_index*3+i]-FBold[i])/massB[IB_index]/2.0 ;
            posB[IB_index*3+i] += items[IDX_dt]*(3.0*velB[IB_index*3+i]-velBold[i])/2.0 ;
        }


        // float ang[3], qua[4] ; 
            // for(int i=0;i<4;i++){qua[i]=quat[IB_index*4+i];}  for(int i=0;i<3;i++){ang[i]=angleVB[IB_index*3+i];}
            // quat[IB_index*4+0] += items[IDX_dt]*(-ang[0]*qua[1] - ang[1]*qua[2] - ang[2]*qua[3])/2.0 ;
            // quat[IB_index*4+1] += items[IDX_dt]*( ang[0]*qua[0] + ang[2]*qua[2] - ang[1]*qua[3])/2.0 ;
            // quat[IB_index*4+2] += items[IDX_dt]*( ang[1]*qua[0] - ang[2]*qua[1] + ang[0]*qua[3])/2.0 ;
            // quat[IB_index*4+3] += items[IDX_dt]*( ang[2]*qua[0] + ang[1]*qua[1] - ang[0]*qua[2])/2.0 ;
        // update angle velocity 
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                angleVB[IB_index*3+i] += items[IDX_dt]
                *(3.0*quaS[IB_index*3+j*3+i]*Torque[IB_index*3+j] - quaSold[j*3+i]*Torqueold[j]) /inertia[IB_index*3+j]/2.0 ; // original
            }
        }        
        // float ang2[3], qua2[4] ; 
            // for(int i=0;i<4;i++){qua2[i]=quat[IB_index*4+i];}  for(int i=0;i<3;i++){ang2[i]=angleVB[IB_index*3+i];}
            // quat[IB_index*4+0] = qua[0] + items[IDX_dt]*(-ang[0]*qua[1] - ang[1]*qua[2] - ang[2]*qua[3] -  ang2[0]*qua2[1] - ang2[1]*qua2[2] - ang2[2]*qua2[3])/4.0  ;
            // quat[IB_index*4+1] = qua[1] + items[IDX_dt]*( ang[0]*qua[0] + ang[2]*qua[2] - ang[1]*qua[3] +  ang2[0]*qua2[0] + ang2[2]*qua2[2] - ang2[1]*qua2[3])/4.0  ;
            // quat[IB_index*4+2] = qua[2] + items[IDX_dt]*( ang[1]*qua[0] - ang[2]*qua[1] + ang[0]*qua[3] +  ang2[1]*qua2[0] - ang2[2]*qua2[1] + ang2[0]*qua2[3])/4.0  ;
            // quat[IB_index*4+3] = qua[3] + items[IDX_dt]*( ang[2]*qua[0] + ang[1]*qua[1] - ang[0]*qua[2] +  ang2[2]*qua2[0] + ang2[1]*qua2[1] - ang2[0]*qua2[2])/4.0  ;


        // update Quaternion & Quaternion
        {float ang[3], qua[4] ; 
            for(int i=0;i<4;i++){qua[i]=quat[IB_index*4+i];}  for(int i=0;i<3;i++){ang[i]=angleVB[IB_index*3+i];}
            quat[IB_index*4+0] += items[IDX_dt]*3.0*(-ang[0]*qua[1] - ang[1]*qua[2] - ang[2]*qua[3])/4.0 ;
            quat[IB_index*4+1] += items[IDX_dt]*3.0*( ang[0]*qua[0] + ang[2]*qua[2] - ang[1]*qua[3])/4.0 ;
            quat[IB_index*4+2] += items[IDX_dt]*3.0*( ang[1]*qua[0] - ang[2]*qua[1] + ang[0]*qua[3])/4.0 ;
            quat[IB_index*4+3] += items[IDX_dt]*3.0*( ang[2]*qua[0] + ang[1]*qua[1] - ang[0]*qua[2])/4.0 ;
            for(int i=0;i<3;i++){ang[i]=angleVBold[i];}
            quat[IB_index*4+0] -= items[IDX_dt]*(-ang[0]*quatold[1] - ang[1]*quatold[2] - ang[2]*quatold[3])/4.0 ;
            quat[IB_index*4+1] -= items[IDX_dt]*( ang[0]*quatold[0] + ang[2]*quatold[2] - ang[1]*quatold[3])/4.0 ;
            quat[IB_index*4+2] -= items[IDX_dt]*( ang[1]*quatold[0] - ang[2]*quatold[1] + ang[0]*quatold[3])/4.0 ;
            quat[IB_index*4+3] -= items[IDX_dt]*( ang[2]*quatold[0] + ang[1]*quatold[1] - ang[0]*quatold[2])/4.0 ;
            for(int i=0;i<4;i++){quatold[IB_index*4+i]=qua[i];qua[i]=quat[IB_index*4+i];}
            // for(int i=0;i<4;i++){quat[IB_index*4+i]/=(powf(qua[0],2)+powf(qua[1],2)+powf(qua[2],2)+powf(qua[3],2));}
        } // */
        set_quaternionS(0,quat[IB_index*4+0],quat[IB_index*4+1],quat[IB_index*4+2],quat[IB_index*4+3],quaS) ;
    }
}

__global__ void update_IBpoint(float *items, int IB_index, float *posB, float *velB, float *angleVB, float *quaS, float *posw, float *oposw, float *nBvec, float *onBvec, float *velw){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
        for(int i=0;i<3;i++){
            posw[id_rho*3+i] = posB[IB_index*3+i] ; velw[id_rho*3+i] = velB[IB_index*3+i] ; nBvec[id_rho*3+i] = 0 ;
            for(int j=0;j<3;j++){
                posw[id_rho*3+i] += quaS[IB_index*9+i*3+j]*oposw[id_rho*3+j] ;
                nBvec[id_rho*3+i]+= quaS[IB_index*9+i*3+j]*onBvec[id_rho*3+j] ;
            }
        }
        velw[id_rho*3+0] += quaS[IB_index*9+0*3+0] * (angleVB[IB_index*3+1]*oposw[id_rho*3+2] - angleVB[IB_index*3+2]*oposw[id_rho*3+1]) 
                          + quaS[IB_index*9+0*3+1] * (angleVB[IB_index*3+2]*oposw[id_rho*3+0] - angleVB[IB_index*3+0]*oposw[id_rho*3+2]) 
                          + quaS[IB_index*9+0*3+2] * (angleVB[IB_index*3+0]*oposw[id_rho*3+1] - angleVB[IB_index*3+1]*oposw[id_rho*3+0]) ;
        velw[id_rho*3+1] += quaS[IB_index*9+1*3+0] * (angleVB[IB_index*3+1]*oposw[id_rho*3+2] - angleVB[IB_index*3+2]*oposw[id_rho*3+1]) 
                          + quaS[IB_index*9+1*3+1] * (angleVB[IB_index*3+2]*oposw[id_rho*3+0] - angleVB[IB_index*3+0]*oposw[id_rho*3+2]) 
                          + quaS[IB_index*9+1*3+2] * (angleVB[IB_index*3+0]*oposw[id_rho*3+1] - angleVB[IB_index*3+1]*oposw[id_rho*3+0]) ;
        velw[id_rho*3+2] += quaS[IB_index*9+2*3+0] * (angleVB[IB_index*3+1]*oposw[id_rho*3+2] - angleVB[IB_index*3+2]*oposw[id_rho*3+1]) 
                          + quaS[IB_index*9+2*3+1] * (angleVB[IB_index*3+2]*oposw[id_rho*3+0] - angleVB[IB_index*3+0]*oposw[id_rho*3+2]) 
                          + quaS[IB_index*9+2*3+2] * (angleVB[IB_index*3+0]*oposw[id_rho*3+1] - angleVB[IB_index*3+1]*oposw[id_rho*3+0]) ;
        }
}

__global__ void search_IBlattice(float *items, int IB_index, int *lattice_id, int *neib, float *posx, float *posy, float *posz, float *posw){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_IBPoints]){
        int id=lattice_id[id_rho], near_id=lattice_id[id_rho] ;
        float dist1 = powf(posx[id]-posw[id_rho*3+0],2) +powf(posy[id]-posw[id_rho*3+1],2) +powf(posz[id]-posw[id_rho*3+2],2) ;
        for(int k=1;k<items[IDX_Q];k++){
            int nei_id= neib[id*(int)items[IDX_Q]+k] ;
            float dist2 = powf(posx[nei_id]-posw[id_rho*3+0],2) +powf(posy[nei_id]-posw[id_rho*3+1],2) +powf(posz[nei_id]-posw[id_rho*3+2],2) ;
            near_id = near_id*(dist1<=dist2) + nei_id*(dist2<dist1) ;
            dist1   = dist1  *(dist1<=dist2) + dist2 *(dist2<dist1) ;
        }
        lattice_id[id_rho] = near_id ;
    }
}


template void IB_csv<float>(int, vector<float>&, vector<float>&, vector<float>&, vector<float>&);
template void IB_csv<double>(int,vector<double>&,vector<double>&,vector<double>&,vector<double>&);
template void set_quaternionS<float> (int IB_index, float q0, float q1, float q2, float q3, vector<float>& quaS) ;
template void set_quaternionS<double>(int IB_index, double q0, double q1, double q2, double q3, vector<double>& quaS) ;
template __global__ void SPM<float>(float*, float, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void SPM<double>(double*, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
template __global__ void SPM_ellipse<float>(float*, float, float, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void get_IBMGw2<float>(float*, int*, int*, float*, float*, float*, float*, float*, float*, float*, float*,float*, float*, float*, float*, float*, float*, float*, float*, float);
template __global__ void get_IBMGw2<double>(double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double);
template __global__ void update_velIBM<float>(float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
template __global__ void update_velIBM<double>(double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
