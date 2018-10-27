#include "Instance.hpp"
#include "lp.hpp"
#include <unordered_set>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
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
    std::vector<std::unordered_set<int>> variables_pack;

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
};
