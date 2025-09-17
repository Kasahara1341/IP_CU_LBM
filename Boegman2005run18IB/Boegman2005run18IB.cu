
#include "../all.hpp"

__device__ float profile_s22(float xi, float r){
    float result =    (0 * ( r > xi/2.0 )) 
                    + (0.5*(sin(3.14159*r/xi) + 1.0) * ( fabsf(r) <= xi/2.0 )) 
                    + (1.0 * ( r < -xi/2.0 )) ;
    return result ;
}
__global__ void SPM_ellipse3D(float *items, float *f, float *posx, float *posy, float *posz, float *velx, float *vely, float *velz, float *velB){
    // smoothed-profile method
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho * (int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        float distance = (0.15*posx[id_rho]+posz[id_rho]-0.29)/sqrtf(1+powf(0.15,2)) ;
        float fx=0, fy=0, fz=0 ;
        fx = (velB[0]*0 - velx[id_rho])/items[IDX_dt] * profile_s22(items[IDX_dz],distance) ;
        fy = (velB[1]*0 - vely[id_rho])/items[IDX_dt] * profile_s22(items[IDX_dz],distance) ;
        fz = (velB[2]*0 - velz[id_rho])/items[IDX_dt] * profile_s22(items[IDX_dz],distance) ;
        for(int k =0;k<items[IDX_Q];k++){
            float tmp = items[IDX_w(k)]*items[IDX_dt] * 3.0
            *( items[IDX_cx(k)]*fx + items[IDX_cy(k)]*fy + items[IDX_cz(k)]*fz )/(powf(items[IDX_c],2)) ;
            f[id_f+k] += tmp ;
            velx[id_rho] += items[IDX_cx(k)] * tmp ;
            vely[id_rho] += items[IDX_cy(k)] * tmp ;
            velz[id_rho] += items[IDX_cz(k)] * tmp ;
        }
    }
}
template<typename Typ>
__global__ void update_velIBM2(Typ *items, Typ *f, Typ *velx, Typ *vely, Typ *velz, Typ *FIBx, Typ *FIBy, Typ *FIBz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    if(id_rho<items[IDX_num_calc]){
        int id_f = id_rho * (int)items[IDX_Q] ;
        for(int k =0;k<items[IDX_Q];k++){
            Typ tmp = items[IDX_dt]*items[IDX_w(k)]*(items[IDX_cx(k)]*FIBx[id_rho] + items[IDX_cy(k)]*FIBy[id_rho] + items[IDX_cz(k)]*FIBz[id_rho]) ;
            f[id_f+k]    += tmp ;
            velx[id_rho] += items[IDX_cx(k)] * tmp ;
            vely[id_rho] += items[IDX_cy(k)] * tmp ;
            velz[id_rho] += items[IDX_cz(k)] * tmp ;
        } 
    }
}
float sech(float x){
    return 2.0/(exp(x)+exp(-x)) ;
}


int main (void){

    // h1 is upper layer. h2 is lower layer
    float H=0.29, h1, h2, eta0, c0, alpha, thicness=0.01 ; 
    h1=H*0.50 ; h2=H-h1 ; eta0 =  h1*0.82 ;
    c0 = sqrt(9.81*20/1020*(h1*h2)/H) ;
    alpha = 1.5 * c0 * (h1-h2)/(h1*h2) ;
    float Ti = 2*6.0/c0, Ts = 6.0/(alpha*eta0) ;
    float slope_eta0 = eta0 / 3.0 , slope= 3.0/20.0 ;

    int i, j, l, k ;
    int *cx, *cy, *cz ;
    vector<int> neib, nextK, nextB, divx, divy ; // devided x, y
    vector<float> f, g, h, Fk, vel_x, vel_y, vel_z, rho, sal, phi, pressure, Fx, Fy, Fz, posx, posy, posz, delX, delY ;
    vector<float> tau, taus ;
    vector<float> item ;
    Items items ; input_items(items,"./input/initial") ; 
    items.PFthick = 3.5*items.dx ; items.sigma = 0.072 *0;
    string Boussinesq_approxi = "on" ; 
    bool Boussi_flag = (strcmp(Boussinesq_approxi.c_str(),"on") ==0) ;
    items.setc(9) ; 
    vector<float> M((int)pow(items.num_velocity,2)), MM((int)pow(items.num_velocity,2)),
    M_inv((int)pow(items.num_velocity,2)), S(items.num_velocity) ;
    set_M<float>(items.num_velocity, M, S, M_inv, MM) ;

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    items.save_interval = 1.0/items.dt ; items.total_count= 75/items.dt ;
    items.save_interval = items.total_count/150 ;
    // items.total_count=20 ; items.save_interval=1 ;
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
        
    cx = new int[items.num_velocity] ; cy= new int[items.num_velocity] ; cz = new int[items.num_velocity] ; 
    for(k=0;k<items.num_velocity;k++){
        cx[k]=items.cx[k]/items.c ; cy[k]=items.cy[k]/items.c ; cz[k]=items.cz[k]/items.c ;
    }

    // divide x, y direction
    for(i=0;i<items.nx;i++){
        float x = (i+0.5)*items.dx*items.ratiox ;
        if( x < 3.0 ){
            if( x < 2.5){
                if(x< 2.0){
                    divx.push_back(30) ; continue ;
                }
                divx.push_back(10) ; continue ;
            }
            divx.push_back(6) ; continue ;
            // divx.push_back(1) ;
        }
        else{divx.push_back(1);}
    } // */
    // divide x, y direction
    // for(i=0;i<items.nx;i++){
        // divx.push_back(1) ;
    // }

    for(j=0;j<items.ny;j++){divy.push_back(1) ; }
    // set wall infomation
    vector<int> lnum ;
    items.num_calc=0 ;
    for(i=0;i<items.nx;i++){
        for(int divi=0;divi<divx[i];divi++){
            for(j=0;j<items.ny;j++){
                for(int divj=0;divj<divy[j];divj++){
                    for(l=0;l<items.nz;l++){
                        float x = ( i + (2.0*divi+1.0)/(2.0*divx[i]) )*items.dx*items.ratiox ;
                        float y = ( j + (2.0*divj+1.0)/(2.0*divy[j]) )*items.dx*items.ratioy ;
                        float z = ( l + 0.5 )*items.dx ;
                        // if(z > 0.29 -6*items.dx - slope * x){
                        if(x > 0){
                            lnum.push_back(items.num_calc) ;
                            items.num_calc+=1 ;
                        }
                        else{lnum.push_back(-1) ; }
                    }
                }
            }
        }
    }    

    printf("set initial condition \n") ;
    float rhoL=1, rhoH=1000.0, muL=1.*pow(10,-5), muH=1.016*pow(10,-3) ;
    { int int_tmp=0 ;
    for(i=0;i<items.nx;i++){
        for(int divi=0;divi<divx[i];divi++){
            for(j=0;j<items.ny;j++){
                for(int divj=0;divj<divy[j];divj++){
                    for(l=0;l<items.nz;l++){
                        if(lnum[int_tmp]<0){ int_tmp+=1 ; continue ;}
                        float local_x = (i + (2*divi+1.0)/(2.0*divx[i]) )*items.dx*items.ratiox ;
                        float local_z = (0.5 + l)*items.dx ;
                        posx.push_back(local_x) ; posz.push_back(local_z) ;
                        posy.push_back( (j + (2*divj+1.0)/(2.0*divy[j]) )*items.dx*items.ratioy ) ;
                        delX.push_back( items.dx*items.ratiox/divx[i]) ; delY.push_back( items.dx*items.ratioy/divy[j]) ;
                        sal.push_back(0) ; phi.push_back(1) ;           
                        sal[sal.size()-1] = 12.13-12.13*tanh((local_z-(h2+(local_x-3)*slope_eta0))/thicness) ; // 6m 
                        // sal[sal.size()-1] = 12.13-12.13*tanh((local_z-(h2+(local_x-3)*slope_eta0*0))/thicness) ; // horizontal ditribution
                        if(posz[posz.size()-1] < 0.29 - slope * posx[posx.size()-1]){
                            sal[sal.size()-1] = 0 ;
                        }
                        // phi[phi.size()-1] = 0.5 - 0.5*tanh( (posz[posz.size()-1] - (d0+eta) )/items.PFthick*2.0 ) ; // 
                        rho.push_back(rhoL+phi[phi.size()-1]*(rhoH-rhoL) + sal[sal.size()-1]*0.824493) ;

                        // pressure.push_back(9.81*(items.dx*items.nz-z) *3.0/(pow(items.c,2))) ; // 静水圧を仮定したp*の分布 
                        pressure.push_back(0) ;
                        vel_x.push_back(0) ; vel_y.push_back(0) ; vel_z.push_back(0) ; 
                        Fx.push_back(0) ; Fy.push_back(0) ; Fz.push_back(0) ;
                        tau.push_back( 1.0*3.0*(muL+phi[phi.size()-1]*(muH-muL))/rho[rho.size()-1]/pow(items.c,2)/items.dt+0.5) ;
                        taus.push_back(1.4*3.0*pow(10,-9)/pow(items.c,2)/items.dt+0.5) ;
                        for(k=0;k<items.num_velocity;k++){
                            float tmp = (vel_x[vel_x.size()-1]*items.cx[k] + vel_y[vel_x.size()-1]*items.cy[k] + vel_z[vel_x.size()-1]*items.cz[k])/pow(items.c,2) ;
                            Fk.push_back(0) ; f.push_back(0) ;
                            g.push_back(items.weight[k]*sal[sal.size()-1]*(1 + 3*tmp + 4.5*tmp*tmp - 1.5*(pow(vel_x[vel_x.size()-1],2)+pow(vel_y[vel_x.size()-1],2)+pow(vel_z[vel_x.size()-1],2))/pow(items.c,2) )) ;
                        }
                        int_tmp+= 1 ;
                    }
                }
            }
        }
    } }
    items.nx=0 ; items.ny=0 ; 
    for(i=0;i<divx.size();i++){items.nx+=divx[i];} for(j=0;j<divy.size();j++){items.ny+=divy[j];}
    for(i=0;i<items.nx;i++){for(j=0;j<items.ny;j++){for(l=0;l<items.nz;l++){
        if(lnum[i*items.ny*items.nz+j*items.nz+l]<0){
            continue ;
        }
        for(k=0;k<items.num_velocity;k++){
            if(i+cx[k]>=0 && i+cx[k]<items.nx && j+cy[k]>=0 && j+cy[k]<items.ny && l+cz[k]>=0 && l+cz[k]<items.nz){
                neib.push_back(lnum[i*items.ny*items.nz+j*items.nz+l + cx[k]*items.ny*items.nz + cy[k]*items.nz + cz[k]]) ;
            }
            else{ neib.push_back(-1) ; }
        }
    }}}

    /* set neighbor wall lattice */
    printf("set neighbor wall lattice \n") ;
    set_neibghor_wall(items,lnum,divx,divy,neib,f,g,Fk,pressure,rho,phi,posx,posy,posz,delX,delY,vel_x,vel_y,vel_z) ;
    cout<<"items nx ="<<items.nx<<" items ny="<<items.ny<<endl<<endl;
    cout<<"number of calculation lattice is "<<items.num_calc<<" wall lattice is "<<rho.size()-items.num_calc<<endl; cout<<""<<endl;
    // hydrostatic_pressure(items,Boussi_flag,neib,pressure,rho,f,posz) ;
    for(i=0;i<items.num_calc;i++){
        for(k=0;k<items.num_velocity;k++){
            float tmp = (vel_x[i]*items.cx[k] + vel_y[i]*items.cy[k] + vel_z[i]*items.cz[k])/pow(items.c,2) ;
            f[i*items.num_velocity+k] = items.weight[k]*(pressure[i]+3.0*tmp+4.5*pow(tmp,2)-1.5*(pow(vel_x[i],2)+pow(vel_y[i],2)+pow(vel_z[i],2))/pow(items.c,2)) ;
            g[i*items.num_velocity+k] = items.weight[k]*phi[i]*(1.0+3.0*tmp+4.5*pow(tmp,2)-1.5*(pow(vel_x[i],2)+pow(vel_y[i],2)+pow(vel_z[i],2))/pow(items.c,2)) ;
        }
    }

    // set IBM points
    items.num_IBMpoints = sqrt(pow(0.29,2)+pow(0.29/0.15,2))/items.dx ;
    vector<float> velB, posB, angleV_B, quaternion, quaS, IB, massB, FB, Torque, densB ;
    vector<int> num_IBMpoints, lattice_id ;
    vector<float> posw, Gw, velw, oposw, onB_vec, nB_vec ;
    // decide IB infomation
    num_IBMpoints.push_back(items.num_IBMpoints) ;
    posB.push_back(0) ; 
    posB.push_back(items.dx*items.ny/2.0) ; 
    posB.push_back(0) ;
    quaternion.push_back(1); 
    for(i=0;i<3;i++){
        quaternion.push_back(0); velB.push_back(0) ; angleV_B.push_back(0) ; Torque.push_back(0) ; FB.push_back(0);
    }
    for(i=0;i<9;i++){quaS.push_back(0);}
    set_quaternionS(0,quaternion[0],quaternion[1],quaternion[2],quaternion[3],quaS) ;
    densB.push_back(1.0*1000) ; massB.push_back(densB[0]) ; // density times area(2D)
    IB.push_back(massB[0]) ; IB.push_back(IB[0]) ; IB.push_back(IB[0]) ;

    cout<<"dens="<<densB[0]<<" massB="<<massB[0]<<endl;

    cout<<"set each IB points"<<endl;
    for(k=0;k<items.num_IBMpoints;k++){
        int near_id=0 ;
        oposw.push_back(20.0/sqrt(409.0)*items.dx*k) ;
        oposw.push_back(0.5*items.dx) ; // y
        oposw.push_back(items.nz*items.dx - 0.29/(items.num_IBMpoints-1.0)*k ) ; 
        // 楕円の法線ベクトルを算出　原点を0とする楕円の法線ベクトル成分は(2x/a^2 , 2y/b^2)
        onB_vec.push_back(3.0/sqrt(409.0)) ;
        onB_vec.push_back(0) ; // y
        onB_vec.push_back(20.0/sqrt(409.0)) ;
        for(i=0;i<3;i++){ 
            Gw.push_back(0) ; velw.push_back(0) ;
            nB_vec.push_back(onB_vec[k*3+i]) ; posw.push_back(oposw[k*3+i]) ;
        }
        float dist1 = 100 ;
        for(i=0;i<items.num_calc;i++){
            float dist2 = sqrt(pow(posx[i]-posw[k*3+0],2) +pow(posy[i]-posw[k*3+1],2) +pow(posz[i]-posw[k*3+2],2) ) ;
            if(dist1>dist2){
                dist1 = dist2 ; near_id = i ;
            }
        }
        lattice_id.push_back(near_id) ;
    }
    cout<<" number of IBM Points = "<<items.num_IBMpoints<< " nz= "<< items.nz<<endl ;

    // index show in spreadsheet
    // https://docs.google.com/spreadsheets/d/1wy2RkS1ECD7LtZCgyQAtm0fKckvEZmfZUeOMNmJwrgk/edit?gid=0#gid=0
    item.push_back(items.dx) ; item.push_back(items.dt) ; item.push_back(items.c) ; 
    item.push_back(items.nx) ; item.push_back(items.ny) ; item.push_back(items.nz) ; item.push_back(items.num_velocity) ; // 0~6
    item.push_back(items.ratiox) ; item.push_back(items.ratioy) ; item.push_back(items.PFthick) ; // 7~9
    // num_calc num_wall
    item.push_back(items.num_calc) ; item.push_back(rho.size()-items.num_calc) ;
    // num_IBM_points IBMdx
    item.push_back(items.num_IBMpoints) ; item.push_back(items.dx/2.0) ;
    item.push_back(items.nu) ; item.push_back(pow(10,-9)) ; item.push_back(items.sigma) ;
    item.push_back(items.tau) ; item.push_back(items.taus) ;
    // wall function用の変数を準備
    vector<int> wall1, wall2, wall3, wall4, wall5, wall6 ;
    set_walln(item, neib, wall1, wall2, wall3, wall4, wall5, wall6) ;
    item.push_back(wall1.size()) ; item.push_back(wall2.size()) ; item.push_back(wall3.size()) ; 
    item.push_back(wall4.size()) ; item.push_back(wall5.size()) ; item.push_back(wall6.size()) ;

    for(i=0;i<items.num_velocity;i++){item.push_back(items.weight[i]) ;}
    for(i=0;i<items.num_velocity;i++){item.push_back(items.cx[i]) ;}
    for(i=0;i<items.num_velocity;i++){item.push_back(items.cy[i]) ;}
    for(i=0;i<items.num_velocity;i++){item.push_back(items.cz[i]) ;}

    if     (items.num_velocity==9 ){set_bound2D(item, items.num_calc, neib, nextK, nextB) ;printf("call set_bound2D\n") ;}
    else if(items.num_velocity==27){set_bound3D(item, items.num_calc, neib, nextK, nextB) ;printf("call set_bound3D\n") ;}

    printf("allocate device memory \n");
    int *d_neib, *d_nextB, *d_nextK ;
    float *d_f, *d_ftmp, *d_fout, *d_feq, *d_g ;
    float *d_posx, *d_posy, *d_posz, *d_delX, *d_delY ;
    float *d_rho, *d_u, *d_v, *d_w, *d_sal, *d_phi, *d_pressure, *d_tau, *d_taus ;
    float *d_phiold, *d_uold, *d_vold, *d_wold ;
    float *d_Fk, *d_Fx, *d_Fy, *d_Fz, *d_FIBx, *d_FIBy, *d_FIBz ;
    float *d_items, *d_M, *d_Minv, *d_S, *d_MM; 
    int *d_wall1, *d_wall2, *d_wall3, *d_wall4, *d_wall5, *d_wall6 ;
    // IBM
    int   *d_lattice_id ; 
    float *d_velB, *d_posB ;
    float *d_angleVB, *d_quaternion, *d_quaS, *d_IB, *d_massB ;
    float *d_FB, *d_Torque, *d_densB, *d_posw, *d_oposw, *d_Gw, *d_velw, *d_nBvec, *d_onBvec ;
    float *d_quatold ;
    cuMallocCopy(&d_neib, neib) ; cuMallocCopy(&d_nextB, nextB) ; cuMallocCopy(&d_nextK,nextK) ;
    cuMallocCopy(&d_M, M)        ; cuMallocCopy(&d_MM, MM) ; 
    cuMallocCopy(&d_Minv, M_inv) ; cuMallocCopy(&d_S, S) ; 
    cuMallocCopy(&d_items,item) ;
    cuMallocCopy(&d_f, f) ;  cuMallocCopy(&d_ftmp, f) ; cuMallocCopy(&d_fout, f) ; cuMallocCopy(&d_feq, Fk) ; cuMallocCopy(&d_Fk, Fk) ;
    cuMallocCopy(&d_g, g) ;  
    cuMallocCopy(&d_pressure,pressure) ; cuMallocCopy(&d_sal,sal) ; cuMallocCopy(&d_phi,phi) ; 
    cuMallocCopy(&d_tau,tau) ; cuMallocCopy(&d_taus,taus) ;    
    cuMallocCopy(&d_Fx,Fx) ; cuMallocCopy(&d_Fy,Fy)   ; cuMallocCopy(&d_Fz,Fz) ;
    cuMallocCopy(&d_FIBx,Fx) ; cuMallocCopy(&d_FIBy,Fy)   ; cuMallocCopy(&d_FIBz,Fz) ;
    cuMallocCopy(&d_rho,rho) ; cuMallocCopy(&d_posx,posx) ; cuMallocCopy(&d_posy,posy) ; cuMallocCopy(&d_posz,posz) ;
    cuMallocCopy(&d_delX, delX) ; cuMallocCopy(&d_delY,delY) ;
    cuMallocCopy(&d_u,vel_x) ; cuMallocCopy(&d_v,vel_y)   ; cuMallocCopy(&d_w,vel_z) ;
    cuMallocCopy(&d_phiold,phi) ; cuMallocCopy(&d_uold,vel_x) ; cuMallocCopy(&d_vold,vel_y) ; cuMallocCopy(&d_wold,vel_z) ;
    cuMallocCopy(&d_wall1, wall1) ; cuMallocCopy(&d_wall2, wall2) ; cuMallocCopy(&d_wall3, wall3) ;
    cuMallocCopy(&d_wall4, wall4) ; cuMallocCopy(&d_wall5, wall5) ; cuMallocCopy(&d_wall6, wall6) ;
    // ibm
    cuMallocCopy(&d_lattice_id,lattice_id) ; cuMallocCopy(&d_velB,velB) ; cuMallocCopy(&d_posB,posB) ; 
    cuMallocCopy(&d_angleVB,angleV_B) ; cuMallocCopy(&d_quaternion,quaternion) ; cuMallocCopy(&d_quaS,quaS) ;
    cuMallocCopy(&d_IB,IB) ; cuMallocCopy(&d_massB,massB) ; cuMallocCopy(&d_FB,FB) ; 
    cuMallocCopy(&d_Torque,Torque) ; cuMallocCopy(&d_densB,densB) ; cuMallocCopy(&d_posw,posw) ; 
    cuMallocCopy(&d_oposw,oposw) ; cuMallocCopy(&d_Gw,Gw) ; cuMallocCopy(&d_velw,velw) ; 
    cuMallocCopy(&d_nBvec,nB_vec) ; cuMallocCopy(&d_onBvec,onB_vec) ; cuMallocCopy(&d_quatold,quaternion) ; 

    //////////////////////////////////////////////////////////////////////////////////////////////////
    output(item,posx,posy,posz,delX,delY,pressure,vel_x,vel_y,vel_z,sal,phi,rho,Fx,Fy,Fz,0,items.save_interval) ;
    IB_csv(0,item, posw, velw, Gw) ;
    printf("start main calculation \n");
    int blockSize = 64;
    int numBlocks = (rho.size() + blockSize - 1) / blockSize ;     
    auto start=chrono::high_resolution_clock::now() ;
    for(int timestep=1 ; timestep<items.total_count+1 ; timestep++){
        // velocity field
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 1, 0, 0, wall1.size(), d_wall1, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 1, 0, 0, wall3.size(), d_wall3, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 0, 1, wall5.size(), d_wall5, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 0, 1, wall6.size(), d_wall6, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 1, 0, wall2.size(), d_wall2, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 1, 0, wall4.size(), d_wall4, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ; // */
        equ_f         <float> <<<numBlocks, blockSize>>>(d_items, d_feq, d_pressure, d_u, d_v, d_w) ;
        Force         <float> <<<numBlocks, blockSize>>>(d_items, Boussi_flag, d_neib, d_f, d_feq, d_tau, d_Fk, d_Fx, d_Fy, d_Fz, d_pressure, d_rho, d_sal, d_phi, d_u, d_v, d_w, d_delX, d_delY, d_posx, d_posy, d_posz) ;
        col_f_MRT     <float> <<<numBlocks, blockSize>>>(d_items, d_tau, d_f, d_ftmp, d_feq, d_Fk, d_M, d_Minv, d_S, d_MM) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_f,d_feq,d_ftmp,d_fout,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ; // 0 => slip ; 1 => bounce back noslip

        // salinity 
        col_g_reg     <float> <<<numBlocks, blockSize>>>(d_items, d_taus, d_g, d_ftmp, d_feq, d_sal, d_u, d_v, d_w) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_g,d_feq,d_ftmp,d_fout,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ;
        // Phase Field
        /*col_PF        <float>       <<<numBlocks, blockSize>>>(d_items, d_neib, d_taus, d_g, d_gtmp, d_geq, d_phi, d_u, d_v, d_w, d_phiold, d_uold, d_vold, d_wold, d_posx, d_posy, d_posz) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_g,d_feq,d_ftmp,d_fout,d_items_adv,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ; // */
        
        update_scalar <float> <<<numBlocks, blockSize>>>(d_items, d_g, d_sal) ;
        // update_scalar <float> <<<numBlocks, blockSize>>>(d_items, d_, d_phi) ;
        update_rho    <float> <<<numBlocks, blockSize>>>(d_items, rhoL, rhoH, d_f, d_Fx, d_Fy, d_Fz, d_pressure, d_sal, d_phi, d_rho, d_u, d_v, d_w) ; 
        LES           <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_tau, d_taus, d_phi, d_rho, muL, muH, d_u, d_v, d_w, d_posx, d_posy, d_posz) ;

        // resetF<float><<<numBlocks, blockSize>>>(d_items, d_Fx, d_Fy, d_Fz, Fx.size()) ;
        for(i=0;i<1;i++){
            // SPM_ellipse3D         <<<numBlocks, blockSize>>>(d_items,d_f,d_posx,d_posy,d_posz,d_u,d_v,d_w,d_velB) ;
            // get_IBMGw2    <float> <<<numBlocks, blockSize>>>(d_items,d_lattice_id,d_neib,d_f,d_tau,d_posx,d_posy,d_posz,d_posw,d_posB,d_nBvec,d_u,d_v,d_w,d_velw,d_Fx,d_Fy,d_Fz,d_Gw,rhoH) ;
            // update_velIBM <float> <<<numBlocks, blockSize>>>(d_items,d_lattice_id,d_f,d_ftmp,d_pressure,d_tau,d_u,d_v,d_w,d_uold,d_vold,d_wold,d_Fx,d_Fy,d_Fz) ;

            get_IBMGw2    <float> <<<numBlocks, blockSize>>>(d_items,d_lattice_id,d_neib,d_f,d_tau,d_posx,d_posy,d_posz,d_posw,d_posB,d_nBvec,d_u,d_v,d_w,d_velw,d_FIBx,d_FIBy,d_FIBz,d_Gw,rhoH) ;
            update_velIBM2<float> <<<numBlocks, blockSize>>>(d_items,d_f,d_u,d_v,d_w,d_FIBx,d_FIBy,d_FIBz) ;
        } //
        set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_rho) ;
        set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_phi) ;
        // set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_u) ;  set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_v) ; set_wall_rho<float><<<numBlocks, blockSize>>>(d_items, d_neib, d_w) ;
        // CUDAのエラーをcheck
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
            return 1;
        } // */
        if(timestep%items.save_interval==0){
            cudaMemcpy(pressure.data(), d_pressure , pressure.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_x.data(), d_u    ,  vel_x.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_y.data(), d_v    ,  vel_y.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_z.data(), d_w    ,  vel_z.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fx.data()   , d_Fx   ,  Fx.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fy.data()   , d_Fy   ,  Fy.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fz.data()   , d_Fz   ,  Fz.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(sal.data()  , d_sal  ,  sal.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(phi.data()  , d_phi  ,  phi.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(rho.data()  , d_rho  ,  rho.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(tau.data()  , d_tau  ,  tau.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            printf("loop%02d  time=%f\n",timestep/items.save_interval,timestep*items.dt) ;
            output<float>(item,posx,posy,posz,delX,delY,pressure,vel_x,vel_y,vel_z,sal,phi,rho,Fx,Fy,Fz,timestep,items.save_interval) ;
            if(isnan(vel_x[0])!=0){
                cout<<"######################################"<<endl<<"Not a number is detected !"
                <<endl<<"######################################"<<endl; break;} // check NAN
            cudaMemcpy(Gw.data()   , d_Gw   ,  Gw.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(posw.data() , d_posw ,  posw.size()  * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(velw.data() , d_velw ,  velw.size()  * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(lattice_id.data(), d_lattice_id , lattice_id.size()* sizeof(int), cudaMemcpyDeviceToHost) ;
            IB_csv(timestep/items.save_interval,item, posw, velw, Gw) ;
        }
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_Fx, d_Fy, d_Fz, Fx.size()) ;
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_FIBx, d_FIBy, d_FIBz, Fx.size()) ;
    }
    cout<<" dz= "<<items.dx<< " dt= " << items.dt<< " nu= "<<items.nu<<" tau= "<<items.tau<<endl;
    cout<<"taus= "<<items.taus<<" ratiox= "<<item[7]<<endl;
    cout<<"nx= "<<items.nx<< " ny= "<<items.ny<< " nz= "<<items.nz<<" num_velocity= "<<items.num_velocity<<endl;
    cout<<"number of calculation lattice is "<<items.num_calc<<" wall lattice is "<<rho.size()-items.num_calc<<endl; cout<<""<<endl;
    auto end=chrono::high_resolution_clock::now() ;
    chrono::duration<float> duration=end-start ;
    cout<<endl;
    cout<<" dz="<<items.dx<<"m, rx="<<items.ratiox<<", ry="<<items.ratioy<<
    ", nx="<<items.nx<<", ny="<<items.ny<<", nz="<<items.nz<<", dt="<<items.dt<<endl;
    /*for(i=50*items.nz;i<50*items.nz+3;i++){
        cout<<i<<" pressure="<<pressure[i]<<endl;
        for(k=0;k<items.num_velocity;k++){
            cout<<k<<" f="<<f[i*items.num_velocity+k]<< "  eq="<<pressure[i]*items.weight[k]<<" Delta="
            <<f[i*items.num_velocity+k]-pressure[i]*items.weight[k]<<" Fk="<<Fk[i*items.num_velocity+k]<< endl;
        }
        cout<<endl;
    } // */
    cout<<"###############################################"<<endl;
    cout<<"compute time = " << duration.count() <<endl;
    cout<<"###############################################"<<endl;

}