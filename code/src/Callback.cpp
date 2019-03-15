#ifndef CALLBACK_HPP
#define CALLBACK_HPP

#include "Flow.hpp"
#include <cmath>
#include <map>
#include <vector>
#include "gurobi_c++.h"

using namespace std;

class Callback: public GRBCallback
{
  public:
    std::vector< std::vector< std::map< int, std::map< int, std::map< int, GRBVar > > > > > vars;
    GRBModel* model;
    int n;
    const Instance &inst_;
    std::vector<std::vector<std::vector<std::vector<GRBVar>>>> enter_flow;
    std::vector<std::vector<std::vector<std::vector<GRBVar>>>> exit_flow;
    std::vector<std::vector<std::vector<std::vector<GRBVar>>>> proccess;
    GRBVar cIdx_;
    Callback(GRBModel *m, const std::vector< std::vector< std::map< int, std::map< int, std::map< int, GRBVar > > > > > &xvars, 
            const Instance &_inst,
            const std::vector<std::vector<std::vector<std::vector<GRBVar>>>> &enter,
            const std::vector<std::vector<std::vector<std::vector<GRBVar>>>> &exit,
            const std::vector<std::vector<std::vector<std::vector<GRBVar>>>> &proc,
            const GRBVar &cIdx
                ) : inst_(_inst) {
      vars = xvars;
      enter_flow = enter;
      exit_flow = exit;
      proccess = proc;
      model = m;
      cIdx_ = cIdx;
    }
  protected:
    void callback() {
      try {
        //cout << "inicio callback. where: " << where << endl;
        if (where == GRB_CB_MIPNODE) { // currently in MIP
            cout << "teste" << endl;
            double objbnd = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
            int bound = ceil(objbnd);
            cout << "Bound: " << bound << endl;
            for (int j = 0; j < inst_.n(); j++){
                // vector< int > idx;
                // vector< double > coef;
                GRBLinExpr expr = 0;
                //idx.push_back( cIdx_ );
                // coef.push_back( 1.0 );
                expr += cIdx_;
                for (int t = 1; t < inst_.maxTime(); t++){
                    
                    for (GRBVar var : enter_flow[j][inst_.m()+1+j][t]){
                        expr -= max(t,bound)*var;
                        // idx.push_back( var );
                        // coef.push_back( -t );
                    }
                }
                // cout << "remove fin("+to_string(j+1)+")"  << endl;
                //model->remove(model->getConstrByName("fin("+to_string(j+1)+")"));
                // model->update();
                //addCut(expr,GRB_GREATER_EQUAL, 0.0);
                //model->addConstr(expr, GRB_GREATER_EQUAL,  0.0, "fin("+to_string(j+1)+")");
                // model->update();
                //model->write("teste_cb.lp");
                // cout << "adicionado fin("+to_string(j+1)+")"  << endl;
                //lp_add_row( mip, idx, coef, "fin("+to_string(j+1)+")", 'G', 0 );
            }
        //   // Found an integer feasible solution - does it visit every node?
        //   double **x = new double*[n];
        //   int *tour = new int[n];
        //   int i, j, len;
        //   for (i = 0; i < n; i++)
        //     x[i] = getSolution(vars[i], n);

        //   if (len < n) {
        //     // Add subtour elimination constraint
        //     GRBLinExpr expr = 0;
        //     for (i = 0; i < len; i++)
        //       for (j = i+1; j < len; j++)
        //         expr += vars[tour[i]][tour[j]];
        //     addLazy(expr <= len-1);
        //   }

        //   for (i = 0; i < n; i++)
        //     delete[] x[i];
        //   delete[] x;
        //   delete[] tour;
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
      
    }
};
#endif