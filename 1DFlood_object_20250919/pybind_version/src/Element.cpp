#include "all.hpp" 

void Element::set_time_solver(shared_ptr<Time_solver> solver){
    time_evo = solver ;
}

void Element::solve_mass_equation(double dt){
    time_evo->update_depth(*this, dt) ;
}

double Element::calc_increment(){
    double fluxin = 0, fluxout = 0, dh ;
    for(auto& node_ptr : up_nodes){
        fluxin  += node_ptr->get_flux() ;
    }
    for(auto& node_ptr : dn_nodes){
        fluxout += node_ptr->get_flux() ;
    }

    dh = (fluxin - fluxout) / (width*length) ;
    return dh ;
}

void Element::set_up_node(shared_ptr<Node> up_node){
    up_nodes.push_back(up_node) ;
}
void Element::set_dn_node(shared_ptr<Node> dn_node){
    dn_nodes.push_back(dn_node) ;
}void Element::set_depth(double setDepth){
    depth = setDepth ;
}

double Element::get_depth()      { return depth     ;}
double Element::get_position()   { return position  ;}
double Element::get_length()     { return length    ;}
double Element::get_elevation()  { return elevation ;}
double Element::get_width()      { return width     ;}
double Element::get_n_manning()  { return n_manning ;}
