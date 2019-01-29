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
#include <unordered_set>

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
}
#endif

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
   // std::vector < std::string > x11_color_;

    // void create_x11();
    void createCompleteGraphDot();

    int manual_cuts();
    double lifting(double c, int *idxs, double *coefs);
    void lifting_linear(int *idxs, double *coefs);
    void lifting_binario(int *idxs, double *coefs);
    int qtd_manual_cuts = 0;
    int qtd_cortes = 0;
    bool clique = true;
    bool continuo = true;
    bool binario = false;


    std::vector< std::vector< std::vector< int > > > xIdx_;
    std::vector< std::vector< std::vector< int > > > eIdx_;
    std::vector< std::vector< int > > fIdx_;

    std::vector<std::vector<std::vector<std::vector<int>>>> process;
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_flow;

    std::vector< int > fim;
    std::vector< std::string > names;

    std::unordered_set<std::vector<int>> variables_pack;


    int getXidx(int j, int m0, int t0, int mf, int tf) const;
    int teto(double v);

};
#endif