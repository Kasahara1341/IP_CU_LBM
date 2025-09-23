#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <type_traits>
#include <iostream>
#include <istream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
using namespace std;

void readCSV(const std::string&, std::vector<double>& , int) ;
void writeCSV(const std::string&, const std::vector<std::vector<double>>&) ;


class Node ; class Time_solver ;

class Element{
    private:
        double position, length, elevation, n_manning, width, depth ;
        vector<shared_ptr<Node>> up_nodes, dn_nodes ;
        unique_ptr<Time_solver> time_evo ;

    public:
    Element(double position_, double length_, double elevation_, double n_manning_, double width_)
        :   position(position_), length(length_), elevation(elevation_), n_manning(n_manning_), width(width_) {}
    void set_time_solver(unique_ptr<Time_solver> solver) ;
    void solve_mass_equation(double) ;
    double calc_increment() ;
    void set_up_node(shared_ptr<Node>) ;
    void set_dn_node(shared_ptr<Node>) ;
    void set_depth(double) ;
    double get_depth() ;
    double get_position() ;
    double get_length() ;
    double get_elevation() ;
    double get_width() ;
    double get_n_manning() ;
} ;

class Node{
    private:
        double flow_q ;
        shared_ptr<Element> up_element, dn_element ;
    public:
    void solve_momentum_equation() ;
    void set_up_element(shared_ptr<Element>) ;
    void set_dn_element(shared_ptr<Element>) ;
    void set_flow_q(double) ;
    shared_ptr<Element> get_up_element() ;
    shared_ptr<Element> get_dn_element() ;
    double get_flux() ;
} ;

class Time_solver{
    public:
        virtual ~Time_solver() = default ;
        // 0は派生クラスが絶対に持たなくてはいけないことを意味する
        virtual void update_depth(Element& , double dt) = 0 ;
} ;

class Euler : public Time_solver{
    public:
        void update_depth(Element& element, double dt) override ;
} ;

class Runge_Kutta_4th : public Time_solver{
    private:
        int stage = 0 ;
        double depth_old, increment[4] ;
    public:
        void update_depth(Element& element, double dt) override ;
        void update_stage() ;
        void set_depth_old(double) ;
        void update_stage0_variables(Element&, double) ;
        void update_stage1_variables(Element&, double) ;
        void update_stage2_variables(Element&, double) ;
        void update_stage3_variables(Element&, double) ;
} ;

class Runge_Kutta_6th : public Time_solver{
    private:
        int stage = 0 ;
        double depth_old, increment[6] ;
    public:
        void update_depth(Element& element, double dt) override ;
        void update_stage() ;
        void set_depth_old(double) ;
        void update_stage0_variables(Element&, double) ;
        void update_stage1_variables(Element&, double) ;
        void update_stage2_variables(Element&, double) ;
        void update_stage3_variables(Element&, double) ;
        void update_stage4_variables(Element&, double) ;
        void update_stage5_variables(Element&, double) ;
} ;