#include "src/all.hpp"
int main(){
    int i, j, number_of_stage=1, output_count=1 ;
    double init_dt, dt, maxt ;
    // inputフォルダから河川の地形情報の読み込み
    vector<double> width1, width2, zb1, zb2, position1, position2, Qin1, Qin2 ;
    readCSV("input/width.csv"    ,width1,0) ;
    readCSV("input/width2.csv"   ,width2,0) ;
    readCSV("input/zb.csv"       ,zb1,0) ;
    readCSV("input/zb2.csv"      ,zb2,0) ;
    readCSV("input/position.csv" ,position1,0) ;
    readCSV("input/position2.csv",position2,0) ;
    readCSV("input/Boundary2.csv",Qin1,1) ;
    readCSV("input/Boundary2.csv",Qin2,2) ;
    init_dt = 0.5 ; 
    dt      = init_dt ;
    maxt    = Qin1.size() ;

    // スマートポインタを用いてオブジェクトvector生成
    vector<shared_ptr<Element>> elements, elements2 ;
    vector<shared_ptr<Node>> nodes, nodes2 ;
    // 本川(高津川)の作成  push_backは引数をvectorの後ろに追加する処理
    for(i=0; i<width1.size();i++){
        elements.push_back(make_shared<Element>(
            position1[i],150,zb1[i],0.03,width1[i])) ;
        nodes.push_back(make_shared<Node>()) ;
    }
    // 支川(匹見川)の作成
    for(i=0; i<width2.size();i++){
        elements2.push_back(make_shared<Element>(
            position2[i],150,zb2[i],0.03,width2[i])) ;
        nodes2.push_back(make_shared<Node>()) ;
    }

    // element & nodeの接続
    // 本川
    for(i=0;i<elements.size();i++){
        nodes[i]   ->set_up_element(elements[i]) ;
        elements[i]->set_dn_node(nodes[i]) ;
    }
    for(i=1;i<elements.size();i++){
        nodes[i-1] ->set_dn_element(elements[i]) ;
        elements[i]->set_up_node(nodes[i-1]) ;
    }
    // 支川
    for(i=0;i<elements2.size();i++){
        nodes2[i]   ->set_up_element(elements2[i]) ;
        elements2[i]->set_dn_node(nodes2[i]) ;
    }
    for(i=1;i<elements2.size();i++){
        nodes2[i-1] ->set_dn_element(elements2[i]) ;
        elements2[i]->set_up_node(nodes2[i-1]) ;
    }
    // 境界条件
    // 上流端境界(Node)
    vector<shared_ptr<Node>> bc_upnode, bc_upnode2 ;
    bc_upnode.push_back( make_shared<Node>()) ; elements[0] ->set_up_node(bc_upnode[0])  ;
    bc_upnode2.push_back(make_shared<Node>()) ; elements2[0]->set_up_node(bc_upnode2[0]) ;
    // 下流端境界(element)
    vector<shared_ptr<Element>> bc_dnelement ;
    bc_dnelement.push_back(make_shared<Element>(
        position1[position1.size()-1]-150, 150,
        zb1[zb1.size()-1]-(zb1[zb1.size()-1]-zb1[0])/(position1[position1.size()-1]-position1[0])*150
        ,0.03,width1[width1.size()-1]
    )) ;
    nodes[nodes.size()-1] -> set_dn_element(bc_dnelement[0]) ;
    // 支川と本川の接続　# 96番目のelementが接続部
    nodes2[nodes2.size()-1] -> set_dn_element(elements[16]) ;
    elements[16]->set_up_node(nodes2[nodes2.size()-1]) ;

    // elementsとnodesをそれぞれ統合(solv_mass_equationなどをまとめたいため 好み)
    vector< vector< shared_ptr<Element>>> combined_elem={elements,elements2} ;
    vector< vector< shared_ptr<Node>>>    combined_node={nodes,nodes2} ;
    // 時間解法の選択
    for(const auto& sublist : combined_elem){
        for(const auto& element_ptr : sublist){
            // Euler
            // element_ptr -> set_time_solver(make_unique<Euler>()) ; number_of_stage=1;
            // RKTablesを指定することで用意した次数のRunge-Kuttaをsetする
            element_ptr->set_time_solver(std::make_unique<Runge_Kutta>(RKTables::RK6())); number_of_stage=6 ;
            // 初期条件
            element_ptr -> set_depth(0) ;
        }
    }
    // 流量初期条件
    for(const auto& sublist : combined_node){
        for(const auto& node_ptr : sublist){
            node_ptr -> set_flow_q(0) ;
        }
    }

    double time=0 ;
    maxt = 40 ;
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////  メイン計算  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
    // 出力用変数 2次元配列のようなもので要素数は[maxt+1][zb1.size()]のようになる
    vector<vector<double>> Q_result(maxt+1,vector<double>(zb1.size())), H_result(maxt+1,vector<double>(zb1.size())), zb_result(maxt+1,vector<double>(zb1.size())) ;
    // 解析時間計測開始
    auto start=chrono::high_resolution_clock::now() ;
    for(int time_increment=1; time/3600<maxt ; time_increment++){
        // 上流端境界条件処理
        bc_upnode[0]    -> set_flow_q(Qin1[(int)(time/3600)]) ;
        bc_upnode2[0]   -> set_flow_q(Qin2[(int)(time/3600)]) ;
        // 下流端境界条件処理
        bc_dnelement[0] -> set_depth(elements[elements.size()-1]->get_depth()) ;
        
        // 水深・流量の計算
        for(int Rkstage=0; Rkstage < number_of_stage ; Rkstage++){
            // #pragma omp parallel for
            // 今回は本川と支川をそれぞれ作ってまとめたので combined_elem.size()=2
            for(int i=0; i < combined_elem.size(); i++) {
                for(int j = 0; j < combined_elem[i].size(); j++) {
                    combined_elem[i][j]->solve_mass_equation(dt);
                }
            }
            // #pragma omp parallel for
            for (int i=0; i < combined_node.size(); i++) {
                for (int j=0; j < combined_node[i].size(); j++) {
                    combined_node[i][j]->solve_momentum_equation();
                }
            }
        }

        // dtの調整
        if(Qin1[(int)time/3600]>200){
            dt = 0.1 ;
        }
        else{
            dt = init_dt ;
        }
        
        // 出力時間を切りよくするための調整
        if(3600*output_count < time+dt){
            dt = 3600*output_count - time ;
            if(dt==0){dt=1;}
        }

        time += dt ;
        if((int)(time)%3600==0){
            cout<<" time = "<<time/3600<<" [hour]"<<" dt:"<<dt<<endl;
            cout<<" Q[50]= "<<elements[50]->get_depth()<<endl;
            cout<<" Qin1 = "<<bc_upnode[0] -> get_flux()<<endl;
            for(i=0;i<zb1.size();i++){
                Q_result [(int)time/3600][i] = nodes[i]   ->get_flux() ;
                H_result [(int)time/3600][i] = elements[i]->get_depth() + elements[i]->get_elevation() ;
                zb_result[(int)time/3600][i] = zb1[i] ;
            }
            output_count+=1 ;
        }
    }
    auto end=chrono::high_resolution_clock::now() ;
    chrono::duration<double> duration=end-start ;
    // \nは改行を意味する
    cout<<"\ncompute time = " << duration.count() << " 秒\n" <<endl;

    // 結果をcsvファイルに出力
    writeCSV("out/Qs.csv",Q_result) ;
    writeCSV("out/Hs.csv",H_result) ;
    writeCSV("out/zb.csv",zb_result) ;
}