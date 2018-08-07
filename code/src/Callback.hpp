#include "Flow_cuts.hpp"



class Callback: public GRBCallback
{
  public:
    GRBVar** vars;
    GRBModel* model;
    int n;
    Callback(GRBModel *m, GRBVar** xvars, int xn) {
      vars = xvars;
      n    = xn;
      model = m;
    }
  protected:
    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          // Found an integer feasible solution - does it visit every node?
          double **x = new double*[n];
          int *tour = new int[n];
          int i, j, len;
          for (i = 0; i < n; i++)
            x[i] = getSolution(vars[i], n);

          if (len < n) {
            // Add subtour elimination constraint
            GRBLinExpr expr = 0;
            for (i = 0; i < len; i++)
              for (j = i+1; j < len; j++)
                expr += vars[tour[i]][tour[j]];
            addLazy(expr <= len-1);
          }

          for (i = 0; i < n; i++)
            delete[] x[i];
          delete[] x;
          delete[] tour;
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
};
