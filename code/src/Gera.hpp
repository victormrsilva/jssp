#ifndef GERA_HPP
#define GERA_HPP

#include "Instance.hpp"
//#include "Callback.hpp"
#include "lp.hpp"
#include "gurobi_c++.h"
#include <map>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

class Gera
{
public:
    Gera( const Instance &_inst );

    void optimize();

    double execute();

    virtual ~Gera();
private:
    const Instance &inst_;

    int cIdx_;
    LinearProgram *mip;
    std::vector< std::vector< std::map< int, std::map< int, std::map< int, int > > > > > xIdx_; // índice dos arcos (y,m0,t0,mf,tf)
   // std::vector < std::string > x11_color_;

    // void create_x11();
    void createCompleteGraphDot();

      // set of entering Geras
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_Gera;
    // set of exiting Geras
    std::vector<std::vector<std::vector<std::vector<int>>>> exit_Gera;
    // proccessing jobs Geras
    std::vector<std::vector<std::vector<std::vector<int>>>> process;

    int getXidx(int j, int m0, int t0, int mf, int tf) const;
    int teto(double v);

};
#endif