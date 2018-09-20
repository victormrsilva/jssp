#ifndef FLOW_HPP
#define FLOW_HPP

#include "Instance.hpp"
//#include "Callback.hpp"
#include "lp.hpp"
#include "gurobi_c++.h"
#include <map>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

class Flow
{
public:
    Flow( const Instance &_inst );

    

    virtual ~Flow();
private:
    const Instance &inst_;

    int cIdx_;
    LinearProgram *mip;
    std::vector< std::vector< std::map< int, std::map< int, std::map< int, int > > > > > xIdx_; // índice dos arcos (y,m0,t0,mf,tf)
    // std::vector < std::string > x11_color_;

    // void create_x11();
    void createCompleteGraphDot();

    template<class T> bool insere_unico(std::vector<T> &vector, T elemento);


    void optimize();

    void cliques(int *idxs,double *coefs);

    std::vector< int > fim;
    std::vector< std::string > names;

    // set of entering flows
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_flow;
    // set of exiting flows
    std::vector<std::vector<std::vector<std::vector<int>>>> exit_flow;
    // proccessing jobs flows
    std::vector<std::vector<std::vector<std::vector<int>>>> process;

    double lifting(int c, int *idxs, double *coefs);
    
    void lifting_linear(int *idxs, double *coefs);
    void lifting_binario(int *idxs, double *coefs);

    int getXidx(int j, int m0, int t0, int mf, int tf) const;
    int teto(double v);

};
#endif