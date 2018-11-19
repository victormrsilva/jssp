#include "Instance.hpp"
#include "lp.hpp"
#include <unordered_set>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

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
}



class Fernando{
public:
    Fernando( const Instance &_inst );


    virtual ~Fernando();
private:
    const Instance &inst_;

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

    int cliques(int *idxs,double *coefs);

    void cgraph_creation();

    int oddHoles();

    
};
