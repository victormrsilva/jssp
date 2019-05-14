#ifndef FLOW_HPP
#define FLOW_HPP

#include "Instance.hpp"
//#include "Callback.hpp"
#include "lp.hpp"
#include "Hash.cpp"
#include "gurobi_c++.h"
#include <map>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

class Flow
{
public:
    Flow( Instance &_inst );



    virtual ~Flow();
private:
    Instance &inst_;

    int cIdx_;
    LinearProgram *mip;
    std::vector< std::vector< std::map< int, std::map< int, std::map< int, int > > > > > xIdx_; // índice dos arcos (y,m0,t0,mf,tf)
    // std::vector < std::string > x11_color_;

    int qtd_cortes = 0;

    // void create_x11();
    void buildProblem();
    void buildCliqueCuts();

    template<class T> bool insere_unico(std::vector<T> &vector, T elemento);

    CGraph *cgraph;

    void optimize();

    int cliques(int *idxs,double *coefs);

    void cgraph_creation();

    std::vector< int > fim;
    std::vector< std::string > names;

    // set of entering flows
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_flow;
    // set of exiting flows
    std::vector<std::vector<std::vector<std::vector<int>>>> exit_flow;
    // proccessing jobs flows
    std::vector<std::vector<std::vector<std::vector<int>>>> process;

    double lifting(double c, int *idxs, double *coefs);

    void lifting_linear(int *idxs, double *coefs);
    void lifting_binario(int *idxs, double *coefs);

    int getXidx(int j, int m0, int t0, int mf, int tf) const;
    double teto(double v);

    bool clique = false;
    bool continuo = true;
    bool binario = true;

    /* melhorar o lst */
    void combinacao(int job, unsigned int tam, std::vector<int> &vec, std::vector<std::vector<int> > &combinacoes);
    void reduz_lst_kondili(int k_max);

    /*corte de fenchel*/
    bool insertVar(std::vector<S> sol, S var);
    bool backtrack(int j, int op, int ti, std::vector<S> sol);
    template <typename T> bool isSubset(std::vector<T> &A, std::vector<T> &B);
    void enumeracao_fenchel(unsigned int r, const std::vector<S> &vars, int index, std::unordered_set<std::vector<S>> &solutions, std::vector<S> solution);
    bool dominancia(std::vector<S> &vec, std::unordered_set<std::vector<S>> &set);
    

};
#endif 