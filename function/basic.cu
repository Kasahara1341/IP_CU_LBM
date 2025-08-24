#include "../all.hpp"

template<typename Typ>
void output(const vector<Typ>& items, const vector<Typ>& posx, const vector<Typ>& posy, const vector<Typ>& posz, const vector<Typ>& delX, const vector<Typ>& delY,  const vector<Typ>& pressure, const vector<Typ>& velx, const vector<Typ>& vely, const vector<Typ>& velz, const vector<Typ>& sal, const vector<Typ>& phi, const vector<Typ>& rho, const vector<Typ>& Fx, const vector<Typ>& Fy, const vector<Typ>& Fz, int counter, int save_interval){
    int i ;
    double total_p, total_s, total_phi, volume ;
    char file_name[128] ;
    FILE *fp ;
    total_p=0 ; total_s=0 ; total_phi=0 ;
    for(i=0;i<items[IDX_num_calc];i++){
        volume = delX[i]*delY[i]/pow(items[IDX_dz],2) ;
        total_p+=pressure[i]*volume ; total_s+=sal[i]*volume ; total_phi+=phi[i]*volume ;
    }

    printf("sum_pressure =%6.9f salinity =%5.10f phi =%5.10f \n",total_p,total_s,total_phi);
    sprintf(file_name,"data/Fakhari_LBM%04d.csv",counter/save_interval) ;
    fp=fopen(file_name,"w") ;
    fprintf(fp,"num_calc,%d\n",(int)(items[IDX_num_calc])) ;
    fprintf(fp,"x       ,y       ,z       ,delx    ,dely    ,pressure,salinity ,phi     ,rho     ,u       ,v       ,w       ,Fx      ,Fy      ,Fz      \n") ;
    for(i=0;i<items[IDX_num_calc];i++){
        // pressure[i]*=1000.0/3.0*pow(items[IDX_c],2);
        fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        posx[i], posy[i], posz[i], delX[i], delY[i],
        pressure[i], sal[i], phi[i], rho[i], velx[i], 
        vely[i], velz[i], Fx[i], Fy[i], Fz[i] 
        ) ;
    }
    fclose(fp) ;

}
template<typename Typ>
void save_status(const vector<Typ>& items, const vector<Typ>& posx, const vector<Typ>& posy, const vector<Typ>& posz, const vector<Typ>& delX, const vector<Typ>& delY,  const vector<Typ>& pressure, const vector<Typ>& velx, const vector<Typ>& vely, const vector<Typ>& velz, const vector<Typ>& sal, const vector<Typ>& phi, const vector<Typ>& rho, const vector<Typ>& f, const vector<Typ>& g, const vector<int>& neib){
    int i ;
    double total_p, total_s, total_phi, volume ;
    char file_name[128] ;
    FILE *fp ;

    sprintf(file_name,"save_status.csv") ;
    fp=fopen(file_name,"w") ;
    fprintf(fp,"num_calc,%d\n",(int)(items[IDX_num_calc])) ;
    for(i=0;i<items[IDX_num_calc];i++){
        fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        posx[i], posy[i], posz[i], delX[i], delY[i],
        pressure[i], sal[i], phi[i], rho[i], velx[i], 
        vely[i], velz[i]) ;
    } fclose(fp) ;
    sprintf(file_name,"save_FGNeib.csv") ; fp=fopen(file_name,"w") ;
    fprintf(fp,"num_calc,%d\n",(int)(items[IDX_num_calc])) ;
    for(i=0;i<items[IDX_num_calc];i++){
        for(int k=0;k<items[IDX_Q];k++){
            fprintf(fp,",%f,%f,%d", f[i*items[IDX_Q]+k], g[i*items[IDX_Q]+k], neib[i*items[IDX_Q]+k]) ;
        }
        fprintf(fp,"\n") ;
    }fclose(fp) ;

}
void input_items(Items &items, const char *fname){
    char filename[128] ; 
    std::snprintf(filename, sizeof(filename), "%s.csv", fname);  // 安全な結合
    printf("input file name= %s\n",filename) ;
    FILE *fp ; fp = fopen(filename,"r") ;
    fscanf(fp, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d\n", &items.nx, &items.ny, &items.nz,
    &items.dx, &items.dt, &items.ratiox, &items.ratioy, &items.gz, &items.nu, &items.save_interval, &items.total_count);
    items.total_count=1+items.save_interval*items.total_count ; items.nu*=pow(10,-6) ;
    items.nx=(int)(items.nx/items.ratiox) ; items.ny=(int)(items.ny/items.ratioy) ;
    if(items.ny<1){items.ny=1; items.ratioy=1; }
    if(items.ny==1){items.num_velocity=9; items.ratioy=1; } else{items.num_velocity=27;} 

}
void Items::setc(int alpha){
    c= dx/dt ; 
    //c = 1.0 ;
    tau=  1.0*3.0*pow(10,-6)    /c/c/dt +0.5 ; 
    taus= 1.4*3.0*pow(10,-alpha)/c/c/dt +0.5 ;
    //tau=nu*3.0/(dt*pow(c,2)) + 0.5 ;
    //taus=pow(10,-9)*3/(dt*pow(c,2)) +0.5 ;
    //taus=pow(10,-5)*3/(dt*pow(c,2)) +0.5 ;
    if(num_velocity==9){
        for(int k=0;k<9;k++){cx[k]=0;cy[k]=0;cz[k]=0;}
        cx[1]=c ; cx[3]=-c ; cx[5]=c ; cx[6]=-c ; cx[7]=-c ; cx[8]=c ;
        cz[2]=c ; cz[4]=-c ; cz[5]=c ; cz[6]=c ; cz[7]=-c ; cz[8]=-c ;
        weight[0]= 4.0/9.0 ;
        for(int k=1;k<5;k++){
            weight[k]   = 1.0/9.0 ; 
            weight[k+4] = 1.0/36.0 ;
        }
    }
    else if(num_velocity==27){
        cx[1]=c ; cx[3]=-c ; cx[7]=c ; cx[8]=c ; cx[9]=-c ; cx[10]=-c ; cx[11]=c ; cx[12]=c ; cx[13]=-c ; cx[14]=-c ;
        cx[19]=c ; cx[20]=c ; cx[21]=c ; cx[22]=c ; cx[23]=-c ; cx[24]=-c ; cx[25]=-c ; cx[26]=-c ;
        cy[2]=c ; cy[4]=-c ; cy[7]=c ; cy[8]=-c ; cy[9]=c ; cy[10]=-c ; cy[15]=c ; cy[16]=c ; cy[17]=-c ; cy[18]=-c ;
        cy[19]=c ; cy[20]=c ; cy[21]=-c ; cy[22]=-c ; cy[23]=c ; cy[24]=c ; cy[25]=-c ; cy[26]=-c ;
        cz[5]=c ; cz[6]=-c ; cz[11]=c ; cz[12]=-c ; cz[13]=c ; cz[14]=-c ; cz[15]=c ; cz[16]=-c ; cz[17]=c ; cz[18]=-c ;
        cz[19]=c ; cz[20]=-c ; cz[21]=c ; cz[22]=-c ; cz[23]=c ; cz[24]=-c ; cz[25]=c ; cz[26]=-c ;
        weight[0]= 8.0/27.0 ;for(int k=1;k<7;k++){weight[k]=2.0/27.0 ;}
        for(int k=7;k<19;k++){weight[k]=1.0/54.0 ;} for(int k=19;k<27;k++){weight[k]=1.0/216.0 ;}
    }
}

template<typename Typ>
void reset_items_dt(vector<Typ> &items, Typ alpha){
    items[IDX_dt] /= alpha ;
    items[IDX_c]  *= alpha ;
    for(int k=0;k<items[IDX_Q];k++){
        items[IDX_cx(k)] *= alpha ;
        items[IDX_cy(k)] *= alpha ;
        items[IDX_cz(k)] *= alpha ;
    }
    
}

__host__ __device__ void set_k_inv(int *k_inv, int num_velocity){
    if(num_velocity==9){
        k_inv[0]=0 ; k_inv[1]=3 ; k_inv[2]=4 ; k_inv[3]=1 ; k_inv[4]=2 ; k_inv[5]=7 ; k_inv[6]=8 ; k_inv[7]=5 ; k_inv[8]=6 ;
    }
    else{
        k_inv[0 ]=0  ; k_inv[1 ]=3  ; k_inv[2 ]=4  ; k_inv[3 ]=1  ; k_inv[4 ]=2  ; k_inv[5 ]=6  ; k_inv[6 ]=5  ; k_inv[7 ]=10 ; k_inv[8 ]=9  ;
        k_inv[9 ]=8  ; k_inv[10]=7  ; k_inv[11]=14 ; k_inv[12]=13 ; k_inv[13]=12 ; k_inv[14]=11 ; k_inv[15]=18 ; k_inv[16]=17 ; k_inv[17]=16 ; 
        k_inv[18]=15 ; k_inv[19]=26 ; k_inv[20]=25 ; k_inv[21]=24 ; k_inv[22]=23 ; k_inv[23]=22 ; k_inv[24]=21 ; k_inv[25]=20 ; k_inv[26]=19 ;
    }
}
__host__ __device__ void set_cpm(int *cp, int *cm, int *cip, int *cim, int num_velocity){
    if(num_velocity==9){
        cp[1]=1  ; cp[3]=3  ; cp[5]=1  ; cp[6]=3  ; cp[7]=3  ; cp[8]=1  ;
        cm[1]=3  ; cm[3]=1  ; cm[5]=3  ; cm[6]=1  ; cm[7]=1  ; cm[8]=3  ;
    }
    else{
        cp[0 ]=0  ; cp[1 ]=1  ; cp[2 ]=2  ; cp[3 ]=3  ; cp[4 ]=4  ; cp[5 ]=0  ; cp[6 ]=0  ; cp[7 ]=7  ; cp[8 ]=8  ;
        cp[9 ]=9  ; cp[10]=10 ; cp[11]=1  ; cp[12]=1  ; cp[13]=3  ; cp[14]=3  ; cp[15]=2  ; cp[16]=2  ; cp[17]=4  ;
        cp[18]=4  ; cp[19]=7  ; cp[20]=7  ; cp[21]=8  ; cp[22]=8  ; cp[23]=9  ; cp[24]=9  ; cp[25]=10 ; cp[26]=10 ;
        cm[0 ]=0  ; cm[1 ]=3  ; cm[2 ]=4  ; cm[3 ]=1  ; cm[4 ]=2  ; cm[5 ]=0  ; cm[6 ]=0  ; cm[7 ]=10 ; cm[8 ]=9  ;
        cm[9 ]=8  ; cm[10]=7  ; cm[11]=3  ; cm[12]=3  ; cm[13]=1  ; cm[14]=1  ; cm[15]=4  ; cm[16]=4  ; cm[17]=2  ;
        cm[18]=2  ; cm[19]=10 ; cm[20]=10 ; cm[21]=9  ; cm[22]=9  ; cm[23]=8  ; cm[24]=8  ; cm[25]=7  ; cm[26]=7  ;

        cip[19]=8  ; cip[20]=8  ; cip[21]=7  ; cip[22]=7  ; cip[23]=10 ; cip[24]=10 ; cip[25]=9  ; cip[26]=9  ;
        cim[19]=9  ; cim[20]=9  ; cim[21]=10 ; cim[22]=10 ; cim[23]=7  ; cim[24]=7  ; cim[25]=8  ; cim[26]=8  ;
    }
}
__host__ __device__ void set_cxypm(int *cxp, int *cxm, int *cyp, int *cym){
    cxp[0 ]=0 ; cxp[1 ]=1 ; cxp[2 ]=0 ; cxp[3 ]=3 ; cxp[4 ]=0 ; cxp[5 ]=0 ; cxp[6 ]=0 ; cxp[7 ]=1 ; cxp[8 ]=1 ;
    cxp[9 ]=3 ; cxp[10]=3 ; cxp[11]=1 ; cxp[12]=1 ; cxp[13]=3 ; cxp[14]=3 ; cxp[15]=0 ; cxp[16]=0 ; cxp[17]=0 ;
    cxp[18]=0 ; cxp[19]=1 ; cxp[20]=1 ; cxp[21]=1 ; cxp[22]=1 ; cxp[23]=3 ; cxp[24]=3 ; cxp[25]=3 ; cxp[26]=3 ;
    cxm[0 ]=0 ; cxm[1 ]=3 ; cxm[2 ]=0 ; cxm[3 ]=1 ; cxm[4 ]=0 ; cxm[5 ]=0 ; cxm[6 ]=0 ; cxm[7 ]=3 ; cxm[8 ]=3 ;
    cxm[9 ]=1 ; cxm[10]=1 ; cxm[11]=3 ; cxm[12]=3 ; cxm[13]=1 ; cxm[14]=1 ; cxm[15]=0 ; cxm[16]=0 ; cxm[17]=0 ;
    cxm[18]=0 ; cxm[19]=3 ; cxm[20]=3 ; cxm[21]=3 ; cxm[22]=3 ; cxm[23]=1 ; cxm[24]=1 ; cxm[25]=1 ; cxm[26]=1 ;

    cyp[0 ]=0 ; cyp[1 ]=0 ; cyp[2 ]=2 ; cyp[3 ]=0 ; cyp[4 ]=4 ; cyp[5 ]=0 ; cyp[6 ]=0 ; cyp[7 ]=2 ; cyp[8 ]=4 ;
    cyp[9 ]=2 ; cyp[10]=4 ; cyp[11]=0 ; cyp[12]=0 ; cyp[13]=0 ; cyp[14]=0 ; cyp[15]=2 ; cyp[16]=2 ; cyp[17]=4 ;
    cyp[18]=4 ; cyp[19]=2 ; cyp[20]=2 ; cyp[21]=4 ; cyp[22]=4 ; cyp[23]=2 ; cyp[24]=2 ; cyp[25]=4 ; cyp[26]=4 ;
    cym[0 ]=0 ; cym[1 ]=0 ; cym[2 ]=4 ; cym[3 ]=0 ; cym[4 ]=2 ; cym[5 ]=0 ; cym[6 ]=0 ; cym[7 ]=4 ; cym[8 ]=2 ;
    cym[9 ]=4 ; cym[10]=2 ; cym[11]=0 ; cym[12]=0 ; cym[13]=0 ; cym[14]=0 ; cym[15]=4 ; cym[16]=4 ; cym[17]=2 ;
    cym[18]=2 ; cym[19]=4 ; cym[20]=4 ; cym[21]=2 ; cym[22]=2 ; cym[23]=4 ; cym[24]=4 ; cym[25]=2 ; cym[26]=2 ;
}



template<typename Typ>
void read_M(int num_velocity, vector<Typ> &M, char *fname){
    int i, j ;
    ifstream infile(fname) ; string line ;
    for(i=0;i<num_velocity && getline(infile, line);i++){
        stringstream ss(line) ; string cell ;
        for(j=0;j<num_velocity && getline(ss,cell, ',');j++){
            stringstream value_stream(cell);
            Typ val;
            value_stream >> val;
            M[i * num_velocity + j] = val;
        }
    }
}
template<typename Typ>
void set_M(int num_velocity, vector<Typ> &M, vector<Typ> &S, vector<Typ> &M_inv, vector<Typ> &MM){
    int i, j, k ;
    char fname[128], fname1[128], fname2[128] ;
    sprintf(fname ,"input/MRT_D3Q%d.csv" ,num_velocity) ;
    sprintf(fname1,"input/S_D3Q%d.csv"   ,num_velocity) ;
    sprintf(fname2,"input/Minv_D3Q%d.csv",num_velocity) ;
    ifstream infile(fname1) ; string line ; string cell ;
    getline(infile, line)   ; stringstream ss(line) ;
    for(i=0;i<num_velocity && getline(ss,cell,',');i++){
        stringstream value_stream(cell) ;
        Typ val ; value_stream >> val ; S[i] = val ;
    }
    read_M(num_velocity,MM,fname) ; read_M(num_velocity,M_inv,fname2) ;
    if     (num_velocity==9 ){S[7]= 0 ; S[8]=S[7] ;}
    else if(num_velocity==27){for(i=5;i<10;i++){S[i]=0 ;}}
    for(i=0;i<num_velocity;i++){
        for(j=0;j<num_velocity;j++){
            M[i*num_velocity+j]=S[i]*MM[i*num_velocity+j] ;
        }
    }
    for(i=0;i<num_velocity;i++){
        for(j=0;j<num_velocity;j++){
            MM[i*num_velocity+j]=0 ;
            for(k=0;k<num_velocity;k++){
                MM[i*num_velocity+j] += M_inv[i*num_velocity+k]*M[k*num_velocity+j] ;
            }
        }
    }
    read_M(num_velocity,M,fname) ;
}

template<typename Typ>
void set_neibghor_wall(Items &items, vector<int> &lnum, vector<int> &divx, vector<int> &divy, vector<int> &neib, vector<Typ> &f, vector<Typ> &g, vector<Typ> &Fk, vector<Typ> &pressure, vector<Typ> &rho, vector<Typ> &phi, vector<Typ> &posx, vector<Typ> &posy, vector<Typ> &posz, vector<Typ> &delX, vector<Typ> &delY, vector<Typ> &velx, vector<Typ> &vely, vector<Typ> &velz){
    int i, j, k , l, int_tmp=0 ;
    int k_inv[27] ; set_k_inv(k_inv,items.num_velocity) ;
    for(i=0;i<items.nx;i++){
        for(j=0;j<items.ny;j++){
            for(l=0;l<items.nz;l++){
                if(lnum[i*items.ny*items.nz+j*items.nz+l]<0){
                    continue ;
                }
                for(int tmp_int1=1;tmp_int1<items.num_velocity;tmp_int1++){
                    if(neib[lnum[i*items.ny*items.nz+j*items.nz+l]*items.num_velocity+tmp_int1]==-1){
                        neib.push_back(lnum[i*items.ny*items.nz+j*items.nz+l]) ;
                        f.push_back(0) ; g.push_back(0) ; Fk.push_back(0) ;
                        for(k=1;k<items.num_velocity;k++){
                            neib.push_back((neib.size()-k)/items.num_velocity) ;
                            f.push_back(0) ; g.push_back(0) ; Fk.push_back(0) ;
                        }
                        pressure.push_back(pressure[lnum[i*items.ny*items.nz+j*items.nz+l]]) ;
                        rho.push_back(rho[lnum[i*items.ny*items.nz+j*items.nz+l]]) ; phi.push_back(phi[lnum[i*items.ny*items.nz+j*items.nz+l]]) ;
                        posx.push_back(posx[lnum[i*items.ny*items.nz+j*items.nz+l]]) ; 
                        posy.push_back(posy[lnum[i*items.ny*items.nz+j*items.nz+l]]) ; posz.push_back(posz[lnum[i*items.ny*items.nz+j*items.nz+l]]) ;
                        delX.push_back(delX[lnum[i*items.ny*items.nz+j*items.nz+l]]) ; delY.push_back(delY[lnum[i*items.ny*items.nz+j*items.nz+l]]) ;
                        velx.push_back(0) ; vely.push_back(0) ; velz.push_back(0) ;
                        for(k=1;k<items.num_velocity;k++){
                            if(neib[lnum[i*items.ny*items.nz+j*items.nz+l]*items.num_velocity+k]==-1){
                                neib[lnum[i*items.ny*items.nz+j*items.nz+l]*items.num_velocity+k] = items.num_calc+int_tmp ;
                                neib[(items.num_calc+int_tmp)*items.num_velocity+k_inv[k]] = lnum[i*items.ny*items.nz+j*items.nz+l] ;
                        }}
                        int_tmp += 1 ;
                        break ;
                    }
                }
    }}} 
}

template<typename Typ>
void hydrostatic_pressure(Items &items, bool Boussi_flag, vector<int> &neib, vector<Typ> &pressure, vector<Typ> &rho, vector<Typ> &f, vector<Typ> &posz, Typ center_z){
    int i, k ;
    int up=2 ; if(items.num_velocity==27){up=5;}
    for(i=items.num_calc-1;i>=0;i--){
        // pressure[i] = pressure[neib[i*items.num_velocity+up]] + 3*9.81 * items.dx/pow(items.c,2) ;
        if(Boussi_flag){pressure[i] = 3*9.81 * (items.nz/2.0*items.dx-posz[i]) /pow(items.c,2)/1000 ;}
        else{ pressure[i] = 3*9.81 * (center_z-posz[i]) /pow(items.c,2) ; }
        for(k=0;k<items.num_velocity;k++){
            f[i*items.num_velocity+k] = pressure[i] * items.weight[k] 
            + 9.81*items.cz[k]*3*items.dt*items.weight[k]/pow(items.c,2) ;
        }
    }
}

template<typename Typ>
__global__ void set_wall_rho(Typ *items, int *neib, Typ *rho){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        rho[id_rho] = rho[neib[id_f]] ;
    }
}
template<typename Typ>
__global__ void reset_pressure(Typ *items, Typ *pressure, Typ *f, Typ *velx, Typ *vely, Typ *velz, Typ alpha){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        pressure[id_rho] /= alpha*alpha ;
        for(int k=0;k<items[IDX_Q];k++){
            Typ tmp = (items[IDX_cx(k)]*velx[id_rho]  + items[IDX_cy(k)]*vely[id_rho] +  items[IDX_cz(k)]*velz[id_rho])/(items[IDX_c]*items[IDX_c]) ;
            f [id_f+k]  = items[IDX_w(k)] * (pressure[id_rho] + 3.0*tmp + 4.5*tmp*tmp-1.5*(velx[id_rho]*velx[id_rho]+vely[id_rho]*vely[id_rho]+velz[id_rho]*velz[id_rho])/(items[IDX_c]*items[IDX_c])) ;
        }
    }
}
template<typename Typ>
__global__ void reset_wall_fin(Typ *items, Typ *fin, Typ *finy){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho>=items[IDX_num_calc] && id_rho<items[IDX_num_calc]+items[IDX_num_wall]){
        int k ;
        for(k=0;k<items[IDX_Q];k++){
            fin[id_f+k] = 0 ; finy[id_f+k] = 0 ;
        }
    }
}
template<typename Typ>
__global__ void LES(Typ *items, int *neib, Typ *tau, Typ *taus, Typ *phi, Typ *rho, Typ muL, Typ muH, Typ *velx, Typ *vely, Typ *velz, Typ *posx, Typ *posy, Typ *posz){
    int id_rho = blockIdx.x * blockDim.x + threadIdx.x ;
    int id_f = id_rho*(int)items[IDX_Q] ;
    if(id_rho<items[IDX_num_calc]){
        int k, up=2, down=2 ; if((int)items[IDX_Q]==27){up=5; down=6;}
        Typ C_csm, Q, E, S, F_cs ;
        Typ nablaxy, nablaxz, nablaxx, nablayx, nablayy, nablayz, nablazx, nablazy, nablazz ;
        Typ array[9] ;
                // LES //
         for(k=0;k<9;k++){
            array[k]=0.0 ;
        }
        S=0.0 ; 
        nablaxx=0.0 ; nablaxy=0.0 ; nablaxz=0.0 ; 
        nablayx=0.0 ; nablayy=0.0 ; nablayz=0.0 ;
        nablazx=0.0 ; nablazy=0.0 ; nablazz=0.0 ;
        //x
        nablaxx+= (velx[neib[id_f+1]]-velx[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        nablaxy+= (vely[neib[id_f+1]]-vely[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        nablaxz+= (velz[neib[id_f+1]]-velz[neib[id_f+3]])/(posx[neib[id_f+1]]-posx[neib[id_f+3]] +items[IDX_dz]/1000) ;
        //y
        if((int)items[IDX_Q]==27){
            nablayx+= (velx[neib[id_f+2]]-velx[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;
            nablayy+= (vely[neib[id_f+2]]-vely[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;
            nablayz+= (velz[neib[id_f+2]]-velz[neib[id_f+4]])/(posy[neib[id_f+2]]-posy[neib[id_f+4]] +items[IDX_dz]/1000) ;
        }
        //z
        nablazx+= (velx[neib[id_f+up]]-velx[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;
        nablazy+= (vely[neib[id_f+up]]-vely[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;
        nablazz+= (velz[neib[id_f+up]]-velz[neib[id_f+down]])/(posz[neib[id_f+up]]-posz[neib[id_f+down]] +items[IDX_dz]/1000) ;

        array[0]=(nablaxx+nablaxx)/2.0 ; array[1]=(nablaxy+nablayx)/2.0 ;
        array[2]=(nablaxz+nablazx)/2.0 ; array[3]=(nablayx+nablaxy)/2.0 ;
        array[4]=(nablayy+nablayy)/2.0 ; array[5]=(nablayz+nablazy)/2.0 ;
        array[6]=(nablazx+nablaxz)/2.0 ; array[7]=(nablazy+nablayz)/2.0 ;
        array[8]=(nablazz+nablazz)/2.0 ;
        Q=-1.0/2.0*(nablaxx*nablaxx+nablaxy*nablayx+nablaxz*nablazx+
                    nablayx*nablaxy+nablayy*nablayy+nablayz*nablazy+
                    nablazx*nablaxz+nablazy*nablayz+nablazz*nablazz) ;
        E=1.0/(2.0*powf(nablaxx+nablaxy+nablaxz+nablayx+nablayy+nablayz+nablazx+nablazy+nablazz,2)+powf(10,-40)) ;
        /* F_cs=Q/E , C'=0.05 */
        F_cs=Q/E ;
        F_cs = fmaxf(-1.0, fminf(Q/E, 1.0)) ;
        C_csm=0.05*powf(abs(F_cs),3.0/2.0) ;
        for(k=0;k<9;k++){
            S+=sqrtf(2.0*powf(array[k],2)) ;
        }
        /* add nu_SGS and update tau */
        // items[10] = powf(10,-6) + C_csm*S*powf(powf(items[IDX_dz],3)*powf(items[7],2),1/3.0) ;
        tau[id_rho]  = 3.0 * ((muL+phi[id_rho]*(muH-muL))/rho[id_rho]+C_csm*S*powf(items[IDX_dz],2)     )/(items[IDX_c]*items[IDX_c]*items[IDX_dt]) + 0.5 ;
        // tau[id_rho]  = 3.0 * ((muL+phi[id_rho]*(muH-muL))/rho[id_rho])/(items[IDX_c]*items[IDX_c]*items[IDX_dt]) + 0.5 ;
        taus[id_rho] = 3.0 * (powf(10,-9)+C_csm*S*powf(items[IDX_dz],2)/1000)/(items[IDX_c]*items[IDX_c]*items[IDX_dt]) + 0.5 ;
    }
}

template void output(const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, const vector<float>&, int, int) ;
template void output(const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, int, int) ;
template void reset_items_dt(vector<float>&, float) ;
template void reset_items_dt(vector<double>&, double) ;
template void read_M(int , vector<float>&, char*) ;
template void read_M(int , vector<double>&, char*) ;
template void set_M(int, vector<float>&, vector<float>&, vector<float>&, vector<float>&) ;
template void set_M(int, vector<double>&, vector<double>&, vector<double>&, vector<double>&) ;
template void set_neibghor_wall(Items&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&) ;
template void set_neibghor_wall(Items&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&) ;
template void hydrostatic_pressure(Items&, bool, vector<int>&, vector<float>&, vector<float> &, vector<float>&, vector<float>&, float) ;
template void hydrostatic_pressure(Items&, bool, vector<int>&, vector<double>&, vector<double> &, vector<double>&, vector<double>&, double) ;
template __global__ void set_wall_rho<float>(float*, int*, float*) ;
template __global__ void set_wall_rho<double>(double*, int*, double*) ;
template __global__ void reset_pressure<double>(double*, double*, double*, double*, double*, double*, double) ;
template __global__ void reset_pressure<float>(float*, float*, float*, float*, float*, float*, float) ;
template __global__ void reset_wall_fin<float>(float*, float*, float*) ;
template __global__ void reset_wall_fin<double>(double*, double*, double*) ;
template __global__ void LES<float>(float*, int*, float*, float*, float*, float*, float, float, float*, float*, float*, float*, float*, float*) ;
template __global__ void LES<double>(double*, int*, double*, double*, double*, double*, double, double, double*, double*, double*, double*, double*, double*) ;
