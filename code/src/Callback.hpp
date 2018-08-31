#include "Flow_cuts.hpp"
#include <cmath>
#include <map>
#include <vector>
#include "gurobi_c++.h"


class Callback: public GRBCallback
{
  public:
    GRBVar** vars;
    GRBModel* model;
    LinearProgram *mip;
    const Instance &inst_;
    const vector<vector<map<int,map<int,map<int,int>>>>> &xIdx_; // índice dos arcos (y,m0,t0,mf,tf)
    const vector<vector<vector<vector<int>>>> &process_;
    int machines,jobs;
    int n;
    Callback(LinearProgram* lp, const Instance &inst,
     vector<vector<map< int,map< int,map< int, int > > > > > &xIdx, 
     vector<vector<vector<vector<int>>>> &process) :
    inst_(inst),mip(lp), xIdx_(xIdx),process_(process){
      
    }
  protected:
    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          // Found an integer feasible solution - does it visit every node?
          /* variáveis em conflito
              1 = as variáveis x(j,m(i),t',m(i+1),t'+d) e x(j,m(i),t',m(i),t'+1) com t' > t, ou seja, as variáveis que ainda podem processar e as de espera nos tempos
              2 = as variáveis x(j,m(i-1), t',m(i),t+d) com t'+d > t
              3 = as variáveis de processamento para a máquina i no tempo t
          */
          int rows = lp_rows(mip);
          CGraph *cgraph = cgraph_create(rows*2); // todos os vértices de menos o c

            vector<int> conflitos;
            //cgraph_add_node_conflicts(cgraph,cIdx_,&conflitos[0],conflitos.size());
            vector<int> indices_conflitos;
            for (int m = 0; m < inst_.m(); m++){
                for (int j = 0; j < inst_.n(); j++){
                    int m0 = inst_.machine(j,m);
                    int mf = (m == inst_.m()-1 ? inst_.m() : inst_.machine(j,m+1));
                    int dur = inst_.time(j,m0); // duration time for machine m0 in job j
                    for (int t = inst_.est(j,m0); t < inst_.lst(j,m0); t++){
                        //cout << j << " " << m0 << " " << t << " " << mf << " " << t+dur << endl;
                        int idx = xIdx_[j][m0+1][t][mf+1][t+dur];
                        // cout << idx << endl;
                        indices_conflitos.emplace_back(idx);
                        conflitos.clear();
                        //cout << "variavel: " <<idx << " " << names[idx] << endl;
                        // caso 1
                        // cout << "caso 1: " << names[idx] << endl;
                        for (int tf = inst_.est(j,m0); tf < inst_.lst(j,m0); tf++){ 
                            if (t == tf) continue;
                            // cout << names[xIdx_[j][m0+1][tf-1][m0+1][tf]] << " " << names[xIdx_[j][m0+1][tf][mf+1][tf+dur]] << " ";
                            // conflitos.emplace_back(xIdx_[j][m0+1][tf-1][m0+1][tf]);
                            conflitos.emplace_back(xIdx_[j][m0+1][tf][mf+1][tf+dur]);
                        }
                        // cout << endl;
                        // caso 2
                        // cout << "caso 2: " << endl;
                        if (m != 0){

                            int m_anterior = inst_.machine(j,m-1);
                            int dur_anterior =  inst_.time(j,m_anterior);
                            if (t-dur_anterior >= inst_.time(j,m_anterior)){
                                for (int tf = t-dur_anterior+1; tf < t; tf++){
                                    // cout << names[xIdx_[j][m_anterior+1][tf][m0+1][tf+dur_anterior]] << " ";
                                    conflitos.emplace_back(xIdx_[j][m_anterior+1][tf][m0+1][tf+dur_anterior]);
                                }
                            }
                        }
                        // caso 3
                        // cout << endl;
                        // cout << "caso 3: " << endl;
                        for (int j_ = 0; j_ < inst_.n(); j_++){
                            for( int var : process_[j_][m0+1][t]){
                                if (var == idx) continue;
                                // cout << names[var] << " ";
                                conflitos.emplace_back(var);
                            }
                        }
                        // cout << endl << endl;

                        // mostra conflitos
                        // cout << "todos os conflitos: " << endl;
                        // cout << names[idx] << " = ";
                        // for (int var : conflitos){
                        //     cout << names[var] << " ";
                        // }
                        // cout << endl;
                        
                        cgraph_add_node_conflicts(cgraph,idx,&conflitos[0],conflitos.size());
                        // idx = xIdx_[j][m0+1][t][m0+1][t+1];
                        // conflitos.clear();
                        // cgraph_add_node_conflicts(cgraph,idx,&conflitos[0],conflitos.size());
                        //cout << "conflitos adicionados no cgraph" << endl;
                        //getchar();
                    }
                }
            }
            cgraph_save(cgraph,"cgraph.txt");
            cout << indices_conflitos.size() << endl;

            double *x = lp_x(mip);
            double *rc = lp_reduced_cost(mip);

            // for (int i = 0; i < names.size(); i++){
            //     cout << names[i] << " " << i << " " << x[i] << " " << rc[i] << endl;
            // }

            // for (int i : indices_conflitos){
            //     cout << i << " ";
            // }

            cout << endl;
            vector<double> x_conflitos = vector<double>(rows*2);
            vector<double> rc_conflitos = vector<double>(rows*2);;
            for (int i; i < rows; i++){
                x_conflitos[i] = x[i];
                x_conflitos[i+rows] = 1-x[i];
                rc_conflitos[i] = rc[i];
                rc_conflitos[i+rows] = (-1)*rc[i];
            }
            delete[] x;
            delete[] rc;
            CliqueSeparation *clique_sep = clq_sep_create(cgraph);
            cout << "clique_sep ok" << endl;
            clq_sep_set_verbose(clique_sep,'T');
            clq_sep_set_rc(clique_sep,&rc_conflitos[0]);//&rc_conflitos[0]);
            cout << "clique_sep_set_rc ok" << endl;
            getchar();
            clq_sep_separate(clique_sep,&x_conflitos[0]);//&x_conflitos[0]);
            cout << "clique_separate ok" << endl;
            getchar();
            const CliqueSet *cliques = clq_sep_get_cliques(clique_sep);
            clq_set_print(cliques);
            getchar();
            int qtd_cliques = clq_set_number_of_cliques(cliques);
            cout << "qtd de cliques" << qtd_cliques << endl;
            for (int i = 0; i < qtd_cliques; i++){
                const IntSet *clq = clq_set_get_clique(cliques,i);
                // cout << "clique " << i << " tamanho " << clq->size << endl;
                // cout << "elementos: ";
                // for (int j = 0; j < clq->size; j++){
                //     cout << clq->elements[j] << " " << names[indices_conflitos[clq->elements[j]]];
                // }
                // cout << endl;
                getchar();
            }
            getchar();
          // double **x = new double*[n];
          // int *tour = new int[n];
          // int i, j, len;
          // for (i = 0; i < n; i++)
          //   x[i] = getSolution(vars[i], n);

          // if (len < n) {
          //   // Add subtour elimination constraint
          //   GRBLinExpr expr = 0;
          //   for (i = 0; i < len; i++)
          //     for (j = i+1; j < len; j++)
          //       expr += vars[tour[i]][tour[j]];
          //   addLazy(expr <= len-1);
          // }

          // for (i = 0; i < n; i++)
          //   delete[] x[i];
          // delete[] x;
          // delete[] tour;
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
};
