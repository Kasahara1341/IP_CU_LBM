#include "all.hpp"

void Euler::update_depth(Element& element, double dt){
    double new_depth ;
    new_depth = element.get_depth() + dt*element.calc_increment() ;
    element.set_depth(new_depth) ;
}

namespace RKTables {
    const ButcherTable& RK2() {
        static const ButcherTable tbl = {
            // int                    次数(ステージ数)
            2,
            // vector<vector<double>> ステージごとでknに乗算する係数
            // k1を計算する際に使う係数はない(k1=dh/dtのため)ので最初は空のvector
            {{       },
            // k2 = depth_old + dt * 1.0*k1
             {1.0    }},
            // vector<double>         最後の更新でknに乗算する係数
            { 0.5, 0.5}
        };
        return tbl;
    }
    const ButcherTable& RK3() {
        static const ButcherTable tbl = {
            3,
            {{                },
             {    0.5         },
             {   -1.0,     2.0}},
            { 1.0/6.0, 2.0/3.0, 1.0/6.0}
        };
        return tbl;
    }
    const ButcherTable& RK4() {
        static const ButcherTable tbl = {
            4,
            {{                  },
             { 0.5              },
             {   0,   0.5       },
             {   0,     0,   1.0}},
            {1.0/6, 1.0/3, 1.0/3, 1.0/6}
        };
        return tbl;
    }
    const ButcherTable& RK6() {
        static const ButcherTable tbl = {
            6,
            {{                                                                                        },
            {           0.2                                                                           },
            {      3.0/40.0,    9.0/40.0,                                                             },
            {           0.3,        -0.9,           1.2                                               },
            {-11.0/54.0,2.5,  -70.0/27.0,     35.0/27.0                                               },
            {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0,  253.0/4096.0              }},
            {    37.0/378.0,           0,   250.0/621.0,      125.0/594.0,             0, 253.0/4096.0}
        };
        return tbl;
    }
}
void Runge_Kutta::update_stage_variables(Element& element, double dt){
    double uppdated_depth ;
    increments[stage] = element.calc_increment() ;
    uppdated_depth = depth_old ;
    for(int i=0;i<stage;i++){
        uppdated_depth += dt * tbl.stage_weights[stage][i] * increments[i] ;
    }
    element.set_depth(uppdated_depth) ;
    // stage発展 or リセット
    update_stage() ;
}
void Runge_Kutta::update_stage(){
    if(stage<tbl.stages-1){ stage += 1 ;}
    else{ stage = 0 ;}
}

void Runge_Kutta::update_depth(Element& element, double dt){
    if(stage==0){   // set old depth
        depth_old = element.get_depth() ;
    }
    update_stage_variables(element,dt) ;
    // 係数が揃ったら次の時間ステップの値に更新する
    if(stage==tbl.stages-1){
        double uppdated_depth=0 ;
        for(int i=0;i< tbl.stages ;i++){
            uppdated_depth += increments[i]*tbl.final_weights[i] ;
        }
        element.set_depth(depth_old + uppdated_depth*dt) ;
    }
}


