#include "all.hpp"

void Node::solve_momentum_equation(){
    double up_H, dn_H, up_depth, dn_depth, up_elev, dn_elev ;
    up_depth = up_element->get_depth()     ; dn_depth = dn_element->get_depth() ;
    up_elev  = up_element->get_elevation() ; dn_elev = dn_element->get_elevation() ;
    up_H  = up_depth + up_elev ; dn_H = dn_depth + dn_elev ;
    double depth, width, n_manning, length, grad ;
    width     = (up_element->get_width()+dn_element->get_width())/2.0 ;
    n_manning = (up_element->get_n_manning() + dn_element->get_n_manning())/2.0 ;
    length    = (up_element->get_length()+dn_element->get_length())/2.0 ;
    if(up_H >= dn_H){
        depth = up_depth ;
        if( up_elev < dn_elev){
            depth = fmax(0.0, up_depth+up_elev-dn_elev) ;
        }
    }
    else{
        depth = dn_depth ;
        if( dn_elev < up_elev){
            depth = fmax(0.0, dn_depth+dn_elev-up_elev) ;
        }
    }
    depth = fmax(0,depth) ;
    grad  = (up_H - dn_H)/length ;
    double sign = (grad>0) - (grad<0) ;
    flux = width/n_manning * pow(depth,5.0/3.0) * sqrt(fabs(grad))*sign ;
}

void Node::set_up_element(shared_ptr<Element> element){
    up_element = element ;
}
void Node::set_dn_element(shared_ptr<Element> element){
    dn_element = element ;
}
void Node::set_flux(double value){
    flux = value ;
}
shared_ptr<Element> Node::get_up_element(){ return up_element ;}
shared_ptr<Element> Node::get_dn_element(){ return dn_element ;}
double  Node::get_flux(){ return flux ;}
