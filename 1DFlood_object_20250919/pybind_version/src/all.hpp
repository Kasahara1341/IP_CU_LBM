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
        shared_ptr<Time_solver> time_evo ;

    public:
    Element(double position_, double length_, double elevation_, double n_manning_, double width_)
        :   position(position_), length(length_), elevation(elevation_), n_manning(n_manning_), width(width_) {}
    void set_time_solver(shared_ptr<Time_solver> solver) ;
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
        double flux ;
        shared_ptr<Element> up_element, dn_element ;
    public:
    void solve_momentum_equation() ;
    void set_up_element(shared_ptr<Element>) ;
    void set_dn_element(shared_ptr<Element>) ;
    void set_flux(double) ;
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

struct ButcherTable {
    int stages;
    vector<vector<double>> stage_weights ;  // stage coefficients
    vector<double> final_weights;                   // output weights
};
namespace RKTables {
    const ButcherTable& RK2();
    const ButcherTable& RK3();
    const ButcherTable& RK4();
    const ButcherTable& RK6();
}
class Runge_Kutta : public Time_solver{
    private:
        ButcherTable tbl ;
        int stage = 0;
        vector<double> increments ;
        double depth_old ;
    public:
        Runge_Kutta(const ButcherTable& Table)
        : tbl(Table), increments(Table.stages, 0.0) {}    
        void update_depth(Element& element, double) override ;
        void update_stage() ;
        void update_stage_variables(Element&, double) ;
} ;

void compute_all(vector<shared_ptr<Element>>& ,vector<shared_ptr<Node>>& ,double , int) ;