
#include "../all.hpp"

__device__ float profile_s22(float limit_lenght, float distance){
    float r = limit_lenght/distance ;
    float result = (0 * (r<0.95)) + (0.5*(sin(3.14159*(r*5 - 5.0)) + 1.0) * (0.95<=r && r<=1.05)) + (1.0 * (1.05<r)) ;
    return result ;
}
__global__ void SPM_ellipse3D(float *items, float Rada, float Radb, float *quaS, float *posB, float *f, float *tau, float *posx, float *posy, float *posz, float *velx, float *vely, float *velz, float *velB){
    // smoothed-profile method
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho * (int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        float X1, Y1, Z1 ;
        X1 = quaS[0]*(posx[id_rho]-posB[0]) + quaS[3]*(posy[id_rho]-posB[1]) + quaS[6]*(posz[id_rho]-posB[2]) ;
        Y1 = quaS[1]*(posx[id_rho]-posB[0]) + quaS[4]*(posy[id_rho]-posB[1]) + quaS[7]*(posz[id_rho]-posB[2]) ;
        Z1 = quaS[2]*(posx[id_rho]-posB[0]) + quaS[5]*(posy[id_rho]-posB[1]) + quaS[8]*(posz[id_rho]-posB[2]) ;
        float distance = powf(X1/Rada,2) + powf(Y1/Rada,2) + powf(Z1/Radb,2), fx, fy, fz ; 
        fx = profile_s22(1.0,distance) * (velB[0] - velx[id_rho])/items[IDX_dt] ;
        fy = profile_s22(1.0,distance) * (velB[1] - vely[id_rho])/items[IDX_dt] ;
        fz = profile_s22(1.0,distance) * (velB[2] - velz[id_rho])/items[IDX_dt] ;
        for(int k =0;k<items[IDX_Q];k++){
            f[id_f+k] += items[IDX_w(k)]*items[IDX_dt] * 3.0
            *( items[IDX_cx(k)]*fx + items[IDX_cy(k)]*fy + items[IDX_cz(k)]*fz )/(powf(items[IDX_c],2)) ;
        }
    }
}
template<typename Typ>
void out_C_D(const vector<Typ>& C_D, const vector<Typ>& C_time){
    string filename = "C_D.csv" ;
    ofstream ofs(filename) ;
    if(!ofs){cout<<"cannot open file "<<filename<<endl; exit(1);}
    ofs<<"C_D,C_time"<<endl ;
    for(int i=0;i<C_D.size();i++){
        ofs<<C_D[i]<<","<<C_time[i]<<endl ;
    }
    ofs.close() ;
} // output C_D

int main (void){

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
    vector<float> vecx_H, vecy_H, vec_velwx, vec_time ;

    items.nu*=1000 ;


    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    items.save_interval = 1.0/items.dt ; items.total_count= 2.05/items.dt ;
    items.save_interval = items.total_count/20 ;
    // items.total_count=2 ; items.save_interval=2 ;
    cout<<"total count= "<<items.total_count<<" save_interval= "<<items.save_interval<<endl;
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
        
    cx = new int[items.num_velocity] ; cy= new int[items.num_velocity] ; cz = new int[items.num_velocity] ; 
    for(k=0;k<items.num_velocity;k++){
        cx[k]=items.cx[k]/items.c ; cy[k]=items.cy[k]/items.c ; cz[k]=items.cz[k]/items.c ;
    }

    // divide x, y direction
    for(i=0;i<items.nx;i++){
        float x = (i+0.5)*items.dx*items.ratiox ;
        // if( 0.2<x && 0.9>x ){    // not uniform
        if( 0.24<x && 0.60>x &&i<0){ // uniform
            if( 0.22<x && 0.8>x){
                divx.push_back(8) ; continue ;
            }
            divx.push_back(2) ; continue ;
        }
        else{divx.push_back(1);}
    } // */

    // divide x, y direction
    for(i=0;i<items.ny;i++){
        float y = (i+0.5)*items.dx*items.ratioy ;
        // if( 0.21<y && 0.63>y ){ // not uniform
        if( 0.24<y && 0.60>y &&i<0){ // uniform
            if( 0.26<y && 0.58>y){
                divy.push_back(8) ; continue ;
            }
            divy.push_back(2) ; continue ;
            // divx.push_back(1) ;
        }
        else{divy.push_back(1);}
    } // */


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
                        if(x > 0){
                        // if(pow(x-0.42,2) + pow(y-0.42,2) > pow(0.03,2) ){
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
                        // sal[sal.size()-1] = 12.13-12.13*tanh((local_z-(h2+(local_x-3)*slope_eta0))/thicness) ; // 6m 
                        // phi[phi.size()-1] = 0.5 - 0.5*tanh( (posz[posz.size()-1] - (d0+eta) )/items.PFthick*2.0 ) ; // 
                        rho.push_back(rhoL+phi[phi.size()-1]*(rhoH-rhoL) + sal[sal.size()-1]*0.824493) ;

                        // pressure.push_back(9.81*(items.dx*items.nz-z) *3.0/(pow(items.c,2))) ; // 静水圧を仮定したp*の分布 
                        pressure.push_back(0) ;
                        vel_x.push_back(0) ; vel_y.push_back(0.0) ; vel_z.push_back(0.0) ; 
                        Fx.push_back(0) ; Fy.push_back(0) ; Fz.push_back(0) ;
                        // tau.push_back( 1.0*3.0*(muL+phi[phi.size()-1]*(muH-muL))/rho[rho.size()-1]/pow(items.c,2)/items.dt+0.5) ;
                        tau.push_back( items.nu*3.0/(pow(items.c,2)*items.dt) + 0.5) ;
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
    cout<<"items nx ="<<items.nx<<" items ny="<<items.ny<<endl<<endl;
    cout<<"lnum="<<lnum.size()<<endl;
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
    // hydrostatic_pressure(items,Boussi_flag,neib,pressure,rho,f,posz) ;
    cout<<"number of calculation lattice is "<<items.num_calc<<" wall lattice is "<<rho.size()-items.num_calc<<endl; cout<<""<<endl;
    for(i=items.num_calc;i<rho.size();i++){Fx.push_back(0);Fy.push_back(0);Fz.push_back(0);}
    for(i=0;i<items.num_calc;i++){
        for(k=0;k<items.num_velocity;k++){
            float tmp = (vel_x[i]*items.cx[k] + vel_y[i]*items.cy[k] + vel_z[i]*items.cz[k])/pow(items.c,2) ;
            f[i*items.num_velocity+k] = items.weight[k]*(pressure[i]+3.0*tmp+4.5*pow(tmp,2)-1.5*(pow(vel_x[i],2)+pow(vel_y[i],2)+pow(vel_z[i],2))/pow(items.c,2)) ;
            g[i*items.num_velocity+k] = items.weight[k]*phi[i]*(1.0+3.0*tmp+4.5*pow(tmp,2)-1.5*(pow(vel_x[i],2)+pow(vel_y[i],2)+pow(vel_z[i],2))/pow(items.c,2)) ;
        }
    }


    // set IBM points //
    vector<float> oposw ;
    float H_axis = items.nz*items.dx , Radius ;
    Radius = 0.01580/2.0 ;

    // read sphere points from csv file
    int number_of_division = 3 ;
    {char sphere_file[100] ;
    sprintf(sphere_file,"sphere_points_n%d.csv",number_of_division) ;
    std::ifstream file(sphere_file);
    if (!file) {std::cerr << "ファイルを開けません\n";return 1;}    
    double spherePoints;
    while (file >> spherePoints) {
        oposw.push_back(spherePoints*Radius) ;
        if (file.peek() == ',') file.ignore(); // カンマを飛ばす
    }    }
    vector<float> velB, posB, angleV_B, quaternion, quaS, IB, massB, FB, Torque, densB ;
    vector<int> num_IBMpoints, lattice_id(oposw.size()/3) ;
    cout<<"oposw.size() = "<<oposw.size()<<endl ;
    vector<float> posw(oposw.size()), Gw(oposw.size()), velw(oposw.size()), onB_vec(oposw.size()), nB_vec(oposw.size()) ;
    items.num_IBMpoints = oposw.size()/3 ;
    // decide IB infomation
    num_IBMpoints.push_back(items.num_IBMpoints) ;
    posB.push_back(5*Radius) ; 
    posB.push_back(0.5*H_axis) ; 
    posB.push_back(0.5*H_axis) ;
    quaternion.push_back(1); 
    for(i=0;i<3;i++){
        quaternion.push_back(0); velB.push_back(0) ; angleV_B.push_back(0) ; Torque.push_back(0) ; FB.push_back(0);
    }
    for(i=0;i<9;i++){quaS.push_back(0);}
    set_quaternionS(0,quaternion[0],quaternion[1],quaternion[2],quaternion[3],quaS) ;
    densB.push_back(1.2*1000) ; massB.push_back(densB[0] * pow(Radius,3) * 4.0/3.0 * 3.141592) ; // density times area(2D)
    IB.push_back(massB[0]*pow(Radius,2) *2.0/5.0 ) ; 
    IB.push_back(IB[0]) ; IB.push_back(IB[0]) ;

    // set each IB points //
    #pragma omp parallel for private(i,j)
    for(k=0;k<items.num_IBMpoints;k++){
        int near_id=0 ;
        // 楕円の法線ベクトルを算出　原点を0とする楕円の法線ベクトル成分は(2x/a^2 , 2y/b^2)
        onB_vec[k*3+0]=(oposw[k*3+0]/Radius) ;
        onB_vec[k*3+1]=(oposw[k*3+1]/Radius) ;
        onB_vec[k*3+2]=(oposw[k*3+2]/Radius) ;
        float norm_nBvec=sqrt(pow(onB_vec[k*3+0],2)+pow(onB_vec[k*3+1],2)+pow(onB_vec[k*3+2],2)) ;
        for(i=0;i<3;i++){ 
            // Gw[k*3+i]=0 ; velw[k*3+i]=0 ;
            float nbvec= 0, pos= 0 ;
            for(j=0;j<3;j++){
                pos   += quaS[i*3+j]*oposw[k*3+j]   ;
                nbvec += quaS[i*3+j]*onB_vec[k*3+j] ;
            }
            nB_vec[k*3+i]=(nbvec/norm_nBvec) ; posw[k*3+i]=(posB[i]+pos) ;
        }
        float dist1 = 100 ;
        for(i=0;i<items.num_calc;i++){
            float dist2 = sqrt(pow(posx[i]-posw[k*3+0],2) +pow(posy[i]-posw[k*3+1],2) +pow(posz[i]-posw[k*3+2],2) ) ;
            if(dist1>dist2){
                dist1 = dist2 ; near_id = i ;
            }
        }
        lattice_id[k]=(near_id) ;
    }

    cout<<" massB = "<<massB[0]<< " radius= "<< Radius <<" S(cm^2)= " << 3.151492*Radius*Radius *10000 <<endl<<
    " num_calc= "<<items.num_calc<<" num_wall= "<<rho.size()-items.num_calc<<" num_f =" <<f.size()<< " num_Fx= "<<Fx.size() <<
    " H = "<<H_axis<<" a = "<<Radius<<" b = "<<Radius<<" \n "<<" num_IBMpoints= "<<items.num_IBMpoints<<
    endl<<endl;
    // printf(" volume = %f m^3 \n",)

    // index show in spreadsheet
    // https://docs.google.com/spreadsheets/d/1wy2RkS1ECD7LtZCgyQAtm0fKckvEZmfZUeOMNmJwrgk/edit?gid=0#gid=0
    item.push_back(items.dx) ; item.push_back(items.dt) ; item.push_back(items.c) ; 
    item.push_back(items.nx) ; item.push_back(items.ny) ; item.push_back(items.nz) ; item.push_back(items.num_velocity) ; // 0~6
    item.push_back(items.ratiox) ; item.push_back(items.ratioy) ; item.push_back(items.PFthick) ; // 7~9
    // num_calc num_wall
    item.push_back(items.num_calc) ; item.push_back(rho.size()-items.num_calc) ;
    // num_IBM_points IBMdx
    item.push_back(items.num_IBMpoints) ; item.push_back(4.0*Radius/sqrt(10.0+2.0*sqrt(5.0))/float(number_of_division)) ;
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

    cout<<"item[num_cul]="<<item[IDX_num_calc]<<"\n"<<" num_vel= "<<item[IDX_Q]<<" item.size()="<< item.size()<<"\n" ;

    printf("allocate device memory \n");
    int *d_neib, *d_nextB, *d_nextK ;
    float *d_f, *d_ftmp, *d_fout, *d_feq, *d_g ;
    float *d_posx, *d_posy, *d_posz, *d_delX, *d_delY ;
    float *d_rho, *d_u, *d_v, *d_w, *d_sal, *d_phi, *d_pressure, *d_tau, *d_taus ;
    float *d_phiold, *d_uold, *d_vold, *d_wold ;
    float *d_Fk, *d_Fx, *d_Fy, *d_Fz ;
    float *d_items, *d_M, *d_Minv, *d_S, *d_MM; 
    // side wall
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
    cuMallocCopy(&d_Fx,Fx)             ; cuMallocCopy(&d_Fy,Fy)   ; cuMallocCopy(&d_Fz,Fz) ;
    cuMallocCopy(&d_rho,rho) ; cuMallocCopy(&d_posx,posx) ; cuMallocCopy(&d_posy,posy) ; cuMallocCopy(&d_posz,posz) ;
    cuMallocCopy(&d_delX, delX) ; cuMallocCopy(&d_delY,delY) ;
    cuMallocCopy(&d_u,vel_x) ; cuMallocCopy(&d_v,vel_y)   ; cuMallocCopy(&d_w,vel_z) ;
    cuMallocCopy(&d_phiold,phi) ; cuMallocCopy(&d_uold,vel_x) ; cuMallocCopy(&d_vold,vel_y) ; cuMallocCopy(&d_wold,vel_z) ;
    cuMallocCopy(&d_wall1, wall1) ; cuMallocCopy(&d_wall2, wall2) ; cuMallocCopy(&d_wall3, wall3) ;
    cuMallocCopy(&d_wall4, wall4) ; cuMallocCopy(&d_wall5, wall5) ; cuMallocCopy(&d_wall6, wall6) ;

    cuMallocCopy(&d_lattice_id,lattice_id) ; cuMallocCopy(&d_velB,velB) ; cuMallocCopy(&d_posB,posB) ; 
    cuMallocCopy(&d_angleVB,angleV_B) ; cuMallocCopy(&d_quaternion,quaternion) ; cuMallocCopy(&d_quaS,quaS) ;
    cuMallocCopy(&d_IB,IB) ; cuMallocCopy(&d_massB,massB) ; cuMallocCopy(&d_FB,FB) ; 
    cuMallocCopy(&d_Torque,Torque) ; cuMallocCopy(&d_densB,densB) ; cuMallocCopy(&d_posw,posw) ; 
    cuMallocCopy(&d_oposw,oposw) ; cuMallocCopy(&d_Gw,Gw) ; cuMallocCopy(&d_velw,velw) ; 
    cuMallocCopy(&d_nBvec,nB_vec) ; cuMallocCopy(&d_onBvec,onB_vec) ; cuMallocCopy(&d_quatold,quaternion) ; 

    //////////////////////////////////////////////////////////////////////////////////////////////////
    output(item,posx,posy,posz,delX,delY,pressure,vel_x,vel_y,vel_z,sal,phi,rho,Fx,Fy,Fz,0,items.save_interval) ;
    IB_csv(0,item, posw, velw, Gw) ;
    vecx_H.push_back(posB[0]/(H_axis)) ; 
    vecy_H.push_back(posB[1]/(H_axis)) ;
    printf("start main calculation \n");
    int blockSize = 64;
    int numBlocks = (rho.size() + blockSize - 1) / blockSize ;     
    auto start=chrono::high_resolution_clock::now() ;
    for(int timestep=1 ; timestep<items.total_count+1 ; timestep++){
        // velocity field
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 1, 0, 0, wall1.size(), d_wall1, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 1, 0, 0, wall3.size(), d_wall3, d_v, d_w, d_u, d_Fy, d_Fz, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 0, 1, wall5.size(), d_wall5, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        // wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 0, 1, wall6.size(), d_wall6, d_u, d_v, d_w, d_Fx, d_Fy, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 1, 0, wall2.size(), d_wall2, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ;
        wall_function <float> <<<numBlocks, blockSize>>>(d_items, d_delX, d_delY, 0, 1, 0, wall4.size(), d_wall4, d_u, d_w, d_v, d_Fx, d_Fz, d_rho) ; 

        equ_f         <float> <<<numBlocks, blockSize>>>(d_items, d_feq, d_pressure, d_u, d_v, d_w) ;
        Force         <float> <<<numBlocks, blockSize>>>(d_items, Boussi_flag, d_neib, d_f, d_feq, d_tau, d_Fk, d_Fx, d_Fy, d_Fz, d_pressure, d_rho, d_sal, d_phi, d_u, d_v, d_w, d_delX, d_delY, d_posx, d_posy, d_posz) ;
        // col_f_MRT     <float> <<<numBlocks, blockSize>>>(d_items, d_tau, d_f, d_ftmp, d_feq, d_Fk, d_M, d_Minv, d_S, d_MM) ;
        col_f_SRT     <float> <<<numBlocks, blockSize>>>(d_items,d_tau,d_f,d_ftmp,d_feq,d_Fk) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_f,d_feq,d_ftmp,d_fout,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ; // 0 => slip ; 1 => bounce back noslip

        // salinity 
        /*col_g_reg     <float> <<<numBlocks, blockSize>>>(d_items, d_taus, d_g, d_ftmp, d_feq, d_sal, d_u, d_v, d_w) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_g,d_feq,d_ftmp,d_fout,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ; //
        // Phase Field
        /*col_PF        <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_taus, d_g, d_ftmp, d_feq, d_phi, d_u, d_v, d_w, d_phiold, d_uold, d_vold, d_wold, d_posx, d_posy, d_posz) ;
        IP_process(d_items,numBlocks,blockSize,d_neib,d_g,d_feq,d_ftmp,d_fout,d_nextB,d_nextK,d_posx,d_posy,d_delX,d_delY,0) ; //
        // set_wall_boundary1<float> <<<numBlocks, blockSize>>>(d_items, items.ny*items.nz, d_wallin , d_g, d_phi, d_u, d_v, d_w) ;
        // set_wall_boundary1<float> <<<numBlocks, blockSize>>>(d_items, items.ny*items.nz, d_wallout, d_g, d_phi, d_u, d_v, d_w) ; */
        
        // update_scalar <float> <<<numBlocks, blockSize>>>(d_items, d_g, d_sal) ;
        // update_scalar <float> <<<numBlocks, blockSize>>>(d_items, d_g, d_phi) ;
        update_rho    <float> <<<numBlocks, blockSize>>>(d_items, rhoL, rhoH, d_f, d_Fx, d_Fy, d_Fz, d_pressure, d_sal, d_phi, d_rho, d_u, d_v, d_w) ; 
        // update_rho    <float> <<<numBlocks, blockSize>>>(d_items, rhoL, rhoH, d_f, d_Fx, d_Fy, d_Fz, d_pressure, d_sal, d_phi, d_rho, d_uold, d_vold, d_wold) ; 
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_Fx, d_Fy, d_Fz, Fx.size()) ;
        // LES           <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_tau, d_taus, d_phi, d_rho, muL, muH, d_u, d_v, d_w, d_posx, d_posy, d_posz) ;
        for(i=0;i<1;i++){
            // SPM           <float> <<<numBlocks, blockSize>>>(d_items, items.dx*items.nx/10 ,d_posB,d_f,d_ftmp,d_tau,d_posx,d_posz,d_Fx,d_Fy,d_Fz,d_u,d_v,d_w,d_velB) ;
            SPM_ellipse3D         <<<numBlocks, blockSize>>>(d_items,Radius,Radius,d_quaS,d_posB,d_f,d_tau,d_posx,d_posy,d_posz,d_u,d_v,d_w,d_velB) ;
            get_IBMGw2    <float> <<<numBlocks, blockSize>>>(d_items,d_lattice_id,d_neib,d_f,d_tau,d_posx,d_posy,d_posz,d_posw,d_posB,d_nBvec,d_u,d_v,d_w,d_velw,d_Fx,d_Fy,d_Fz,d_Gw) ;
            update_velIBM <float> <<<numBlocks, blockSize>>>(d_items,d_lattice_id,d_f,d_ftmp,d_pressure,d_tau,d_u,d_v,d_w,d_uold,d_vold,d_wold,d_Fx,d_Fy,d_Fz) ;

            update_IBbody    <<<numBlocks, blockSize>>>(d_items,0,d_massB,d_densB,d_IB,d_FB,d_posB,d_Torque,d_velB,d_quaternion,d_quaS,d_angleVB,d_posw,d_Gw,d_quatold) ;
            update_IBpoint   <<<numBlocks, blockSize>>>(d_items,0,d_posB,d_velB,d_angleVB,d_quaS,d_posw,d_oposw,d_nBvec,d_onBvec,d_velw) ;
            search_IBlattice <<<numBlocks, blockSize>>>(d_items,0,d_lattice_id,d_neib,d_posx,d_posy,d_posz,d_posw) ;
        } //
        // set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_rho) ; 
        // set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_phi) ;
        // set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_u) ;  set_wall_rho  <float> <<<numBlocks, blockSize>>>(d_items, d_neib, d_v) ; set_wall_rho<float><<<numBlocks, blockSize>>>(d_items, d_neib, d_w) ;

        
        // CUDAのエラーをcheck
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
            return 1;
        } 
        if(timestep%items.save_interval==0){
            cudaMemcpy(pressure.data(), d_pressure , pressure.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_x.data(), d_u    ,  vel_x.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_y.data(), d_v    ,  vel_y.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(vel_z.data(), d_w    ,  vel_z.size() * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fx.data()   , d_Fx   ,  Fx.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fy.data()   , d_Fy   ,  Fy.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Fz.data()   , d_Fz   ,  Fz.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(rho.data()  , d_rho  ,  rho.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(tau.data()  , d_tau  ,  tau.size()   * sizeof(float),         cudaMemcpyDeviceToHost) ;
            printf("loop%02d  time=%f\n",timestep/items.save_interval,timestep*items.dt) ;
            output<float>(item,posx,posy,posz,delX,delY,pressure,vel_x,vel_y,vel_z,sal,phi,rho,Fx,Fy,Fz,timestep,items.save_interval) ;
            if(isnan(vel_x[items.ny*items.nz*5])!=0){
                cout<<"######################################"<<endl<<"Not a number is detected !"
                <<endl<<"######################################"<<endl; break;} // check NAN
            cudaMemcpy(Gw.data()   , d_Gw   ,  Gw.size()    * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(posw.data() , d_posw ,  posw.size()  * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(velw.data() , d_velw ,  velw.size()  * sizeof(float),         cudaMemcpyDeviceToHost) ;
            cudaMemcpy(lattice_id.data(), d_lattice_id , lattice_id.size()* sizeof(int), cudaMemcpyDeviceToHost) ;
            IB_csv(timestep/items.save_interval,item, posw, velw, Gw) ;
            cudaMemcpy(posB.data(), d_posB , posB.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(velB.data(), d_velB , velB.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(Torque.data(), d_Torque , Torque.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cudaMemcpy(FB.data(), d_FB , FB.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cout<<" position    velocity    Torque  Force"<<endl;
            for(i=0;i<3;i++){
                printf("%f  %f  %f  %f\n",posB[i],velB[i],Torque[i],FB[i]);
            }
            vec_velwx.push_back(velB[0]) ; vec_time.push_back(timestep*items.dt) ;
            float x_H=posB[0]/(H_axis), y_H=posB[1]/(H_axis) ; 
            cout<<" x/H= "<<x_H<<" y/H= "<<y_H<<endl;
            vecx_H.push_back(x_H) ; vecy_H.push_back(y_H) ;

            cudaMemcpy(quaternion.data(), d_quaternion , quaternion.size()* sizeof(float), cudaMemcpyDeviceToHost) ;
            cout<<"Q0= "<<quaternion[0]<<"  Q2= "<< quaternion[2]<<" thetax "<<acos(quaternion[0])*180/3.141592 *2<<" thetaz "<<asin(quaternion[2])*180/3.141592 *2<<   endl;
        }
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_Fx, d_Fy, d_Fz, Fx.size()) ;
        resetF<float><<<numBlocks, blockSize>>>(d_items, d_Gw, d_Gw, d_Gw, Gw.size()) ; 
    }

    out_C_D(vec_velwx,vec_time) ;
    cout<<" dz= "<<items.dx<< " dt= " << items.dt<< " nu= "<<items.nu<<" tau= "<<tau[0]<<endl;
    cout<<"taus= "<<items.taus<<" ratiox= "<<item[7]<<endl;
    cout<<"nx= "<<items.nx<< " ny= "<<items.ny<< " nz= "<<items.nz<<" num_velocity= "<<items.num_velocity<<endl;
    cout<<"number of calculation lattice is "<<items.num_calc<<" wall lattice is "<<rho.size()-items.num_calc<<endl; cout<<""<<endl;
    auto end=chrono::high_resolution_clock::now() ;
    chrono::duration<float> duration=end-start ;
    cout<<endl;
    cout<<" dz="<<items.dx<<"m, rx="<<items.ratiox<<", ry="<<items.ratioy<<
    ", nx="<<items.nx<<", ny="<<items.ny<<", nz="<<items.nz<<", dt="<<items.dt<<endl;

    cout<<"###############################################"<<endl;
    cout<<"compute time = " << duration.count() <<endl;
    cout<<"###############################################"<<endl;
    cout<<"Np = "<<items.num_IBMpoints<<" dx = "<<items.dx<<" IDX_dIBM = "<<item[IDX_dIBM]<< " dIBM/dx = "<<pow(item[IDX_dIBM]/items.dx,1) <<endl; 


}