#include "Instance.hpp"
#include "lp.hpp"

extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

class Fernando{
public:
    Fernando( const Instance &_inst );

    void optimize();

    void cliques(int *idxs,double *coefs);

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

    double lifting(int c, int *idxs, double *coefs);
    void lifting_linear(int *idxs, double *coefs);
    void lifting_binario(int *idxs, double *coefs);

    int cIdx_;
    LinearProgram *mip;
    int teto(double v);
};
