#include "all.hpp"

void Euler::update_depth(Element& element, double dt){
    double new_depth ;
    new_depth = element.get_depth() + dt*element.calc_increment() ;
    element.set_depth(new_depth) ;
}

void Runge_Kutta_4th::update_stage(){
    if(stage<3){ stage += 1 ;}
    else{ stage = 0 ;}
}

void Runge_Kutta_4th::set_depth_old(double value){
    depth_old = value ;
}
//////////////////////////////////////////////////////////////////////////
// RKについての一般化？検討 with chat GPT
void Runge_Kutta_4th::update_stage_variables(Element& element, double dt, int order, int stage){
    double uppdated_depth ;
    if(order == 4){
        if(stage==0){
            double coeff[stage] ;
        }
    }
    increment[0] = element.calc_increment() ;
    uppdated_depth = element.get_depth() + 0.5*increment[0]*dt ;
    element.set_depth(uppdated_depth) ;
    update_stage() ;
    
}
void update_stage(Element& element, int i, double t, double dt) {
    double tmp_depth = element.get_depth();
    // y + sum_{j<i} a_ij * k_j * dt を構成
    for (int j = 0; j < i; j++) {
        tmp_depth += dt * tbl->a[i][j] * increments[j];
    }
    // ステージ i の f(...) を評価
    Element tmp_elem = element;   // コピーして一時的に状態を作る
    tmp_elem.set_depth(tmp_depth);
    increments[i] = tmp_elem.calc_increment();
}
namespace RKTables {
    inline const ButcherTableau& RK2() {
        static const ButcherTableau tbl = {
            2,
            {{   0,   0},
             { 1.0,   0}},
            {  0.5, 0.5}
        };
        return tbl;
    }
    inline const ButcherTableau& RK3() {
        static const ButcherTableau tbl = {
            3,
            {{      0,       0,       0},
             {    0.5,       0,       0},
             {   -1.0,     2.0,       0}},
            { 1.0/6.0, 2.0/3.0, 1.0/6.0}
        };
        return tbl;
    }
    inline const ButcherTableau& RK4() {
        static const ButcherTableau tbl = {
            4,
            {{   0,   0,  0,0},
             { 0.5,   0,  0,0},
             {   0, 0.5,  0,0},
             {   0,   0,1.0,0}},
            {1.0/6, 1.0/3, 1.0/3, 1.0/6}
        };
        return tbl;
    }
    inline const ButcherTableau& RK6() {
        static const ButcherTableau tbl = {
            6,
            {
            {           0.2,           0,             0,                0,            0,             0},
            {      3.0/40.0,    9.0/40.0,             0,                0,            0,             0},
            {           0.3,        -0.9,           1.2,                0,            0,             0},
            {-11.0/54.0,2.5,  -70.0/27.0,     35.0/27.0,                0,            0,             0},
            {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0,             0},
            },
            {    37.0/378.0,           0,   250.0/621.0,      125.0/594.0,             0, 253.0/4096.0}
        };
        return tbl;
    }
}
// */
//////////////////////////////////////////////////////////////////////////

void Runge_Kutta_4th::update_stage0_variables(Element& element, double dt){
    double uppdated_depth ;
    increment[0] = element.calc_increment() ;
    uppdated_depth = element.get_depth() + 0.5*increment[0]*dt ;
    element.set_depth(uppdated_depth) ;
    update_stage() ;
    
}
void Runge_Kutta_4th::update_stage1_variables(Element& element, double dt){
    double uppdated_depth ;
    increment[1] = element.calc_increment() ;
    uppdated_depth = element.get_depth() + 0.5*increment[1]*dt ;
    element.set_depth(uppdated_depth) ;
    update_stage() ;
    
}
void Runge_Kutta_4th::update_stage2_variables(Element& element, double dt){
    double uppdated_depth ;
    increment[2] = element.calc_increment() ;
    uppdated_depth = element.get_depth() + increment[2]*dt ;
    element.set_depth(uppdated_depth) ;
    update_stage() ;
    
}
void Runge_Kutta_4th::update_stage3_variables(Element& element, double dt){
    increment[3] = element.calc_increment() ;
    update_stage() ;
}

void Runge_Kutta_4th::update_depth(Element& element, double dt){
    if(stage==0){
        depth_old = element.get_depth() ;
        update_stage0_variables(element,dt) ;
    }
    else if(stage==1){ update_stage1_variables(element,dt); }
    else if(stage==2){ update_stage2_variables(element,dt); }
    else if(stage==3){ 
        update_stage3_variables(element,dt) ;
        double uppdated_depth=0, coeff[4]={1,2,2,1} ;
        for(int i=0;i<4;i++){
            uppdated_depth += increment[i]*coeff[i] / 6.0 ;
        }
        element.set_depth(depth_old + uppdated_depth*dt) ;
    }
    else{cout<<"Runge Kutta error stage access!!!\n";}
}