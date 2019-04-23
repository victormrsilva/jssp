#ifndef FLOWTESTES_H
#define FLOWTESTES_H

#include "Instance.hpp"
#include "lp.hpp"
#include <unordered_set>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

#ifndef HASH_H
#define HASH_H
namespace std {
    template <>
    struct hash<std::vector<int>> {
        size_t operator()(const vector<int>& v) const {
        std::hash<int> hasher;
        std::size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
        }
    };
        struct hash<std::vector<Flow_testes::S>> {
        size_t operator()(const vector<Flow_testes::S>& v) const {
        std::hash<Flow_testes::S> hasher;
        std::size_t seed = 0;
        for (S i : v) {
            seed ^= hasher(i.var) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
        }
    };
}
#endif


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
    void combinacao(int tam, std::vector<int> &vec, std::vector<std::vector<int> > &combinacoes);

    int cliques(int *idxs,double *coefs);

    void cgraph_creation();

    int oddHoles();

    void makespanProblem();

    struct S {
        int i;
        int j;
        int t;
        int var;
        S& operator =(const S& a)
        {
            i = a.i;
            j = a.j;
            t = a.t;
            var = a.var;
            return *this;
        }

        inline bool operator==(S a) {
            if (a.var==var)
                return true;
            else
                return false;
        }
    };

    std::vector< std::vector<Flow_testes::S> > solutions;

    int maxOperationsBT = 5;

    bool insertVar(std::vector<Flow_testes::S> sol, Flow_testes::S var);
    bool backtrack(int j, int op, int ti, std::vector<Flow_testes::S> sol);
    
    void fenchel(int ti, int tf); // pega as variáveis por um intervalo de tempo para fazer o corte
    
};

#endif