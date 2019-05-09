#ifndef FLOWTESTES_H
#define FLOWTESTES_H

#include "Instance.hpp"
#include "lp.hpp"
#include "Hash.cpp"
#include <unordered_set>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

class Flow_testes{
public:
    Flow_testes( Instance &_inst );
    void elimina_variavel_flow(int k_max);
    void elimina_variavel_compact(int k_max);
    void elimina_variavel_kondili(int k_max);
    void inicioBT();

    virtual ~Flow_testes();


private:
    Instance &inst_;

    std::vector< int > fim;
    std::vector< std::string > names;
    

    std::vector< std::vector< std::vector< int > > > xIdx_;
    std::vector< std::vector< std::vector< int > > > eIdx_;
    std::vector< std::vector< int > > fIdx_;

    std::vector<std::vector<std::vector<std::vector<int>>>> process;
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_flow;

    double lifting(double c, int *idxs, double *coefs);
    void lifting_linear(int *idxs, double *coefs);
    void lifting_binario(int *idxs, double *coefs);
    int manual_cuts();
    int qtd_manual_cuts = 0;
    //std::vector<std::vector<int>> variables_pack;
    std::unordered_set<std::vector<int>> variables_pack;

    int qtd_cortes = 0;
    int cIdx_;
    LinearProgram *mip;
    double teto(double v);
    CGraph *cgraph;
    bool clique = true;
    bool continuo = true;
    bool binario = false;
    void optimize();
    void combinacao(unsigned int tam, std::vector<int> &vec, std::vector<std::vector<int> > &combinacoes);

    int cliques(int *idxs,double *coefs);

    void cgraph_creation();

    int oddHoles();

    void makespanProblem();



    std::vector< std::vector<S> > solutions;

    int maxOperationsBT = 5;

    bool insertVar(std::vector<S> sol, S var);
    bool backtrack(int j, int op, int ti, std::vector<S> sol);
    
    int fenchel(int ti, int tf); // pega as variáveis por um intervalo de tempo para fazer o corte

    template <typename T> bool isSubset(std::vector<T> &A, std::vector<T> &B);

    void enumeracao_fenchel(unsigned int r, const std::vector<S> &vars, int index, std::unordered_set<std::vector<S>> &solutions, std::vector<S> solution);

    bool dominancia(std::vector<S> &vec, std::unordered_set<std::vector<S>> &set);
    
};



#endif