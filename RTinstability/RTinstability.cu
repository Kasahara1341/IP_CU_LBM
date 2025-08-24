
#include "../all.hpp"


int main (void){

    // h1 is upper layer. h2 is lower layer
    float H=0.29, h1, h2, eta0, c0, alpha, thicness=0.005 ; 
    h2=H*0.30 ; h1=H-h2 ; eta0 =  h2*0.45 ;
    c0 = sqrt(9.81*20/1020*(h1*h2)/H) ;
    alpha = 1.5 * c0 * (h1-h2)/(h1*h2) ;
    float Ti = 2*6.0/c0, Ts = 6.0/(alpha*eta0) ;
    float slope_eta0 = eta0 / 3.0 , slope= 3.0/20.0 ;

    int i, j, l, k ;
    int *cx, *cy, *cz ;
    vector<int> neib, nextK, nextB ;
    vector<float> f, g, h, Fk, vel_x, vel_y, vel_z, rho, sal, phi, pressure, Fx, Fy, Fz, posx, posy, posz ;
    vector<float> tau, taus ;
    vector<float> item, items_adv(15) ;
    Items items ; 
    input_items(items,"./input/initial") ; items.PFthick = 5*items.dx ;
    items.setc(9) ; set_Lagrange_coeff<float>(items_adv, items.ratiox, items.ratioy) ;
    vector<float> M((int)pow(items.num_velocity,2)), MM((int)pow(items.num_velocity,2)),
    M_inv((int)pow(items.num_velocity,2)), S(items.num_velocity) ;
    set_M<float>(items.num_velocity, M, S, M_inv, MM) ;

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    items.save_interval = 0.05/items.dt ; items.total_count= 1.5/items.dt ;

    // float tmp=0.01 ; 
    // items.save_interval=1 ; items.total_count=10 ;
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
        
    cx = new int[items.num_velocity] ; cy= new int[items.num_velocity] ; cz = new int[items.num_velocity] ; 
    for(k=0;k<items.num_velocity;k++){
        cx[k]=items.cx[k]/items.c ; cy[k]=items.cy[k]/items.c ; cz[k]=items.cz[k]/items.c ;
    }

    // set wall infomation
    int *lnum=new int[items.nx*items.ny*items.nz] ; // lattice number index
    items.num_calc=0 ;
    for(i=0;i<items.nx;i++){
        for(j=0;j<items.ny;j++){
            for(l=0;l<items.nz;l++){
                float x=(0.5+i)*items.dx*items.ratiox ;
                float y=(0.5+j)*items.dx*items.ratioy ;
                float z=(0.5+l)*items.dx ;
                if(z > 0){
                    lnum[i*items.ny*items.nz+j*items.nz+l]=items.num_calc ;
                    items.num_calc+=1 ;
                }
                else{
                    lnum[i*items.ny*items.nz+j*items.nz+l]=-1 ;
                }
            }
        }
    }    

    printf("set initial condition \n") ;
    float rhoL=300, rhoH=1000.0 ;
    for(i=0;i<items.nx;i++){
        for(j=0;j<items.ny;j++){
            for(l=0;l<items.nz;l++){
                if(lnum[i*items.ny*items.nz+j*items.nz+l]<0){
                    continue ;
                }
                vel_x.push_back(0) ; vel_y.push_back(0) ; vel_z.push_back(0) ; 
                posx.push_back((0.5+i)*items.dx*items.ratiox) ; posy.push_back((0.5+j)*items.dx*items.ratioy) ; posz.push_back((0.5+l)*items.dx) ;
                sal.push_back(0) ; phi.push_back(0) ;           
                // sal[sal.size()-1] = 12.13-12.13*tanh((posz[posz.size()-1]-(h2+(posx[posx.size()-1]-3)*slope_eta0))/thicness) ; // 6m 
                phi[phi.size()-1] = 0.5 + 0.5*tanh(
                    (posz[posz.size()-1]- (items.dx*items.nz/2.0 + 0.1*cos(2*3.1415*posx[posx.size()-1]/(items.dx*items.nx*items.ratiox)) )  )/items.PFthick*2
                ) ; // 6m
                rho.push_back(rhoL+phi[phi.size()-1]*(rhoH-rhoL) + sal[sal.size()-1]*0.824493) ;

                pressure.push_back(1000*9.81*(items.dx*items.nz-posz[posz.size()-1]) *3.0/(1000.0*pow(items.c,2)) ) ; // 静水圧を仮定したp*の分布 
                // pressure.push_back(0) ;

                Fx.push_back(0) ; Fy.push_back(0) ; Fz.push_back(-9.81*1000) ;
                tau.push_back( 1.0*3.0*pow(10,-2)/pow(items.c,2)/items.dt+0.5) ;
                taus.push_back(1.4*3.0*pow(10,-9)/pow(items.c,2)/items.dt+0.5) ;
                for(k=0;k<items.num_velocity;k++){
                    Fk.push_back(3*items.weight[k]*items.dt*(Fz[Fz.size()-1]*items.cz[k])/(items.c*1000)) ;
                    f.push_back(pressure[pressure.size()-1]*items.weight[k]) ;// - Fk[Fk.size()-1]/2.0) ; // f0 = p*wi - Fk
                    Fk[Fk.size()-1] = 0 ;
                    // g.push_back(sal[sal.size()-1]*items.weight[k]) ;
                    g.push_back(phi[phi.size()-1]*items.weight[k]) ;
                    if(i+cx[k]>=0 && i+cx[k]<items.nx && j+cy[k]>=0 && j+cy[k]<items.ny && l+cz[k]>=0 && l+cz[k]<items.nz){
                        neib.push_back(lnum[i*items.ny*items.nz+j*items.nz+l + cx[k]*items.ny*items.nz + cy[k]*items.nz + cz[k]]) ;
                    }
                    else{ neib.push_back(-1) ; }
                }
                Fz[Fz.size()-1] = 0 ;
            }
        }
    } 

    /* set neighbor wall lattice */
    printf("set neighbor wall lattice \n") ;
    set_neibghor_wall(items,lnum,neib,f,g,Fk,rho,phi,posx,posy,posz,vel_x,vel_y,vel_z) ;

    // index show in spreadsheet
    // https://docs.google.com/spreadsheets/d/1wy2RkS1ECD7LtZCgyQAtm0fKckvEZmfZUeOMNmJwrgk/edit?gid=0#gid=0
    item.push_back(items.dx) ; item.push_back(items.dt) ; item.push_back(items.c) ; 
    item.push_back(items.nx) ; item.push_back(items.ny) ; item.push_back(items.nz) ; item.push_back(items.num_velocity) ; // 0~6
    item.push_back(items.ratiox) ; item.push_back(items.ratioy) ; item.push_back(items.PFthick) ; // 7~9
    item.push_back(items.num_calc) ; item.push_back(rho.size()-items.num_calc) ;
    item.push_back(items.nu) ; item.push_back(pow(10,-9)) ;
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
    int *d_neib ;
    int *d_nextB, *d_nextK ;
    float *d_f, *d_ftmp, *d_fout, *d_feq ;
    float *d_g, *d_gtmp, *d_gout, *d_geq ;
    float *d_posx, *d_posy, *d_posz, *d_rho, *d_u, *d_v, *d_w, *d_sal, *d_phi, *d_pressure, *d_tau, *d_taus ;
    float *d_phiold, *d_uold, *d_vold, *d_wold ;
    float *d_Fk, *d_Fx, *d_Fy, *d_Fz ;
    float *d_items, *d_items_adv, *d_M, *d_Minv, *d_S, *d_MM; 
    int *d_wall1, *d_wall2, *d_wall3, *d_wall4, *d_wall5, *d_wall6 ;
    cuMallocCopy(&d_neib, neib) ; cuMallocCopy(&d_nextB, nextB) ; cuMallocCopy(&d_nextK,nextK) ;
    cuMallocCopy(&d_M, M)        ; cuMallocCopy(&d_MM, MM) ; 
    cuMallocCopy(&d_Minv, M_inv) ; cuMallocCopy(&d_S, S) ; 
    cuMallocCopy(&d_items,item) ; cuMallocCopy(&d_items_adv,items_adv) ; 
    cuMallocCopy(&d_f, f) ;  cuMallocCopy(&d_ftmp, f) ; cuMallocCopy(&d_fout, f) ; cuMallocCopy(&d_feq, Fk) ; cuMallocCopy(&d_Fk, Fk) ;
    cuMallocCopy(&d_g, g) ;  cuMallocCopy(&d_gtmp, g) ; cuMallocCopy(&d_gout, g) ; cuMallocCopy(&d_geq, Fk) ;
    cuMallocCopy(&d_pressure,pressure) ; cuMallocCopy(&d_sal,sal) ; cuMallocCopy(&d_phi,phi) ; cuMallocCopy(&d_phiold,phi) ;
    cuMallocCopy(&d_tau,tau) ; cuMallocCopy(&d_taus,taus) ;    
    cuMallocCopy(&d_Fx,Fx)             ; cuMallocCopy(&d_Fy,Fy)   ; cuMallocCopy(&d_Fz,Fz) ;
    cuMallocCopy(&d_rho,rho) ; cuMallocCopy(&d_posx,posx) ; cuMallocCopy(&d_posy,posy) ; cuMallocCopy(&d_posz,posz) ;
    cuMallocCopy(&d_u,vel_x) ; cuMallocCopy(&d_v,vel_y)   ; cuMallocCopy(&d_w,vel_z) ;
    cuMallocCopy(&d_uold,vel_x) ; cuMallocCopy(&d_vold,vel_y) ; cuMallocCopy(&d_wold,vel_z) ;
    cuMallocCopy(&d_wall1, wall1) ; cuMallocCopy(&d_wall2, wall2) ; cuMallocCopy(&d_wall3, wall3) ;
    cuMallocCopy(&d_wall4, wall4) ; cuMallocCopy(&d_wall5, wall5) ; cuMallocCopy(&d_wall6, wall6) ;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    output(item,posx,posy,posz,pressure,vel_x,vel_y,vel_z,phi,Fx,Fy,Fz,0,items.save_interval) ;
    printf("start main calculation \n");
    int blockSize = 64;
    // int numBlocks = (items.num_calc+num_wall + blockSize - 1) / blockSize ;     
    int numBlocks = (rho.size() + blockSize - 1) / blockSize ;     
    auto start=chrono::high_resolution_clock::now() ;
    float dt_ratio = 0.2 ; int dt_change = 1 ;
    for(int timestep=1 ; timestep<items.total_count+1 ; timestep++){
    // for(int timestep=1;timestep<1+2;timestep++){
        if(timestep==1){
            reset_items_dt(item, dt_ratio) ; items.save_interval /= dt_ratio ; reset_pressure<float><<<numBlocks, blockSize>>>(d_items, d_pressure, dt_ratio) ;
            cudaMemcpy(d_items, item.data(), item.size() * sizeof(float), cudaMemcpyHostToDevice) ;
        }
        if(timestep==dt_change*items.save_interval+1){
            reset_items_dt(item, float(1.0/dt_ratio)) ; items.save_interval /= 1.0/dt_ratio ; timestep = dt_change*items.save_interval+1 ;
            reset_pressure<float><<<numBlocks, blockSize>>>(d_items, d_pressure, 1.0/dt_ratio) ;
            cudaMemcpy(d_items, item.data(), item.size() * sizeof(float), cudaMemcpyHostToDevice) ;
        } // */
        /*wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratioy, 1.0         , IDX_wall1_size, d_wall1, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratioy, 1.0         , IDX_wall3_size, d_wall3, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratiox, items.ratioy, IDX_wall5_size, d_wall5, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratiox, items.ratioy, IDX_wall6_size, d_wall6, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratiox, 1.0         , IDX_wall2_size, d_wall2, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, items.ratiox, 1.0         , IDX_wall4_size, d_wall4, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ; // */
        equ_f         <float> <<<numBlocks, blockSize>>>(d_items, d_feq, d_pressure, d_u, d_v, d_w) ;
        Force         <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_Fk, d_Fx, d_Fy, d_Fz, d_pressure, d_rho, d_sal, d_phi, d_posx, d_posy, d_posz) ;
        col_f_MRT     <float> <<<numBlocks, blockSize>>>(d_items, d_tau, d_f, d_ftmp, d_feq, d_Fk, d_M, d_Minv, d_S, d_MM) ;
        // equ_g         <float> <<<numBlocks, blockSize>>>(d_items, d_geq, d_sal, d_u, d_v, d_w) ;
        // col_g_reg     <float> <<<numBlocks, blockSize>>>(d_items, d_taus, d_g, d_gtmp, d_geq) ;
        col_PF        <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_taus, d_g, d_gtmp, d_geq, d_phi, d_u, d_v, d_w, d_phiold, d_uold, d_vold, d_wold, d_posx, d_posy, d_posz) ;
        /*set_wall_f    <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_ftmp) ;
        set_wall_f    <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_gtmp) ;
        set_out       <float> <<<numBlocks, blockSize>>>(d_items, d_items_adv, d_neib, d_ftmp, d_fout, d_feq, d_f) ;
        set_out       <float> <<<numBlocks, blockSize>>>(d_items, d_items_adv, d_neib, d_gtmp, d_gout, d_geq, d_g) ;
        set_fin       <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextK, d_feq, d_f) ;
        set_fin       <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextK, d_geq, d_g) ;
        reset_f       <float> <<<numBlocks, blockSize>>>(d_items, d_ftmp, d_feq, d_f) ;
        reset_f       <float> <<<numBlocks, blockSize>>>(d_items, d_gtmp, d_geq, d_g) ;
        set_wall_f_inv<float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_f) ;
        set_wall_f_inv<float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_fout) ;
        set_wall_f_inv<float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_g) ;
        set_wall_f_inv<float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_nextB, d_nextK, d_gout) ;
        set_tmp       <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_f, d_ftmp, d_fout, d_feq) ;
        set_tmp       <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_g, d_gtmp, d_gout, d_geq) ; // */
        // 0 => slip ; 1 => bounce back noslip,
        set_wall_propagation<float> <<<numBlocks, blockSize>>>(d_items, 0, d_neib, d_nextK, d_nextB, d_ftmp) ;
        set_wall_propagation<float> <<<numBlocks, blockSize>>>(d_items, 0, d_neib, d_nextK, d_nextB, d_gtmp) ;
        propagation   <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_f, d_ftmp) ;
        propagation   <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_g, d_gtmp) ;
        update_scalar <float> <<<numBlocks, blockSize>>>(d_items, d_g, d_phi) ;
        update_rho    <float> <<<numBlocks, blockSize>>>(d_items, rhoL, rhoH, d_f, d_Fx, d_Fy, d_Fz, d_pressure, d_sal, d_phi, d_rho, d_u, d_v, d_w) ; 
        // LES           <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_tau, d_taus, d_u, d_v, d_w, d_posx, d_posy, d_posz) ;
        set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_rho) ; set_wall_rho<float><<<numBlocks, blockSize>>>(d_items, d_neib, d_phi) ;
        reset_wall_fin<float> <<<numBlocks, blockSize>>>(d_items, d_feq, d_f) ; reset_wall_fin<float><<<numBlocks, blockSize>>>(d_items, d_geq, d_g) ;
        // CUDAのエラーをチェ�?ク
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
            return 1;
        } // */
        if(timestep%items.save_interval==0){
            cudaMemcpy(pressure.data()  , d_pressure  ,  items.num_calc * sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_x.data(), d_u    ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_y.data(), d_v    ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_z.data(), d_w    ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fx.data()   , d_Fx   ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fy.data()   , d_Fy   ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fz.data()   , d_Fz   ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            cudaMemcpy(sal.data()  , d_phi  ,  items.num_calc * sizeof(float),           cudaMemcpyDeviceToHost) ;
            printf("loop%02d  time=%f\n",timestep/items.save_interval,timestep*items.dt) ;
            output<float>(item,posx,posy,posz,pressure,vel_x,vel_y,vel_z,sal,Fx,Fy,Fz,timestep,items.save_interval) ;
            cudaMemcpy(g.data()  , d_geq, g.size()   * sizeof(float),cudaMemcpyDeviceToHost) ;
            // cudaMemcpy(rho.data(), d_rho , rho.size() * sizeof(float),cudaMemcpyDeviceToHost) ;
            i=47500 ; j=48080 ; 
            for(k=0;k<items.num_velocity;k++){
                // cout<<"i="<<i<<" g["<<k<<"]="<<g[i*items.num_velocity+k]<<" j="<<j<<" g="<<g[j*items.num_velocity+k]<< endl;
                // cout<<"i="<<j<<" g["<<k<<"]="<<g[j*items.num_velocity+k]<< endl;
                // cout<<neib[i*items.num_velocity+k]<<endl;
            } // */
        }
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_Fx, d_Fy, d_Fz) ;
    }
    cout<<" dz= "<<items.dx<< " dt= " << items.dt<< " nu= "<<items.nu<<" tau= "<<items.tau<<endl;
    cout<<"taus= "<<items.taus<<" ratiox= "<<item[7]<<endl;
    cout<<"nx= "<<items.nx<< " ny= "<<items.ny<< " nz= "<<items.nz<<" num_velocity= "<<items.num_velocity<<endl;
    cout<<"number of calculation lattice is "<<items.num_calc<<" wall lattice is "<<rho.size()-items.num_calc<<endl; cout<<""<<endl;
    cout<<"h1= " <<h1<<" h2= " <<h2<< " eta0= " <<eta0<<endl;
    cout<<"c0= " <<c0<< " alpha= " <<alpha<<endl;
    cout<<"Ti= " <<Ti<< " Ts= " <<Ts<<endl;
    auto end=chrono::high_resolution_clock::now() ;
    chrono::duration<float> duration=end-start ;
    cout<<endl;
    cout<<" dz="<<items.dx<<"m, rx="<<items.ratiox<<", ry="<<items.ratioy<<
    ", nx="<<items.nx<<", ny="<<items.ny<<", nz="<<items.nz<<", dt="<<items.dt<<endl;
    cout<<"compute time = " << duration.count() <<endl<<endl;

    for(i=0;i<wall4.size();i++){
        // cout<<i<<" wall4 "<<wall4[i]<<endl;
    }

}