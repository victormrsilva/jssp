#include "Kondili.hpp"

#include <vector>
#include <string>
#include <cfloat>
#include <iostream>
using namespace std;

Kondili::Kondili( const Instance &_inst ) : inst_(_inst) { // já inicializa a variável privada _inst com o valor passado por referencia
    mip = lp_create();

    // variáveis de decisão
    xIdx_ = vector<vector<vector<int>>>(inst_.n(),vector<vector<int>>(inst_.m(),vector<int>(inst_.maxTime())));
    

    vector< string > names; // nome das variáveis
    vector< double > lb; // lower bound
    vector< double > ub; // upper bound
    vector< double > obj; // se é objetivo?
    vector< char > integer; // variável inteira?

    // criação das variáveis x
    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            for (int t = 0; t < inst_.maxTime(); t++){
                xIdx_[i][j][t] = names.size(); // número do índice da variável (vai de 1 até n*m)
                names.push_back("x("+to_string(i+1)+","+to_string(j+1)+","+to_string(t+1)+")"); // nome dessa variável
                lb.push_back(0.0);
                ub.push_back(DBL_MAX);
                obj.push_back(0.0);
                integer.push_back(1);
            }
        }
    }

    // c var
    cIdx_ = names.size();
    names.push_back("C");
    lb.push_back( 0.0 );
    ub.push_back( DBL_MAX );
    obj.push_back( 1.0 );
    integer.push_back( 1 );

    // adiciona colunas ao solver
    lp_add_cols( mip, obj, lb, ub, integer, names );
    cout << "Number of variables: " << names.size() << endl;
    // restrição de tempo
    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            vector< int > idx;
            vector< double > coef;
            for (int t = 0; t < inst_.maxTime(); t++){

                idx.push_back( xIdx_[i][j][t] );
                coef.push_back( 1.0 );
                // adiciona restrição.
            }
            lp_add_row( mip, idx, coef, "time("+to_string(i+1)+","+to_string(j+1)+")", 'E', 1 );
        }
    }
    cout << "time restriction added" << endl;
   // makespan
    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            vector< int > idx;
            vector< double > coef;
            idx.push_back( cIdx_ );
            coef.push_back( 1 );

            for (int t = 0; t < inst_.maxTime(); t++){
                idx.push_back( xIdx_[i][j][t] );
                coef.push_back( -1 * (t+inst_.time(j,i)) );
                // adiciona restrição.
            }
            lp_add_row( mip, idx, coef, "makespan("+to_string(i+1)+","+to_string(j+1)+")", 'G', 0 );
        }
    }
    cout << "makespan restriction added" << endl;
    // restrição process
    for (int i = 0; i < inst_.m(); i++){
        for (int t = 0; t < inst_.maxTime(); t++){
            vector< int > idx;
            vector< double > coef;

            for (int j = 0; j < inst_.n(); j++){
                //bool mudou = false;
                for (int t_aux = 0; t_aux < inst_.time(j,i); t_aux++){
                    if (t - t_aux < 0){
                        break;
                    }
                    idx.push_back( xIdx_[i][j][t - t_aux] );
                    coef.push_back( 1 );
                }
            }

            lp_add_row( mip, idx, coef, "processing("+to_string(i+1)+","+to_string(t+1)+")", 'L', 1 );
        }
    }
    cout << "processing restriction added" << endl;
    // restrição ord

    for (int j = 0; j < inst_.n(); j++){
        for (int h = 1; h < inst_.m(); h++){

            vector< int > idx;
            vector< double > coef;

            for (int t = 0; t < inst_.maxTime(); t++){
                idx.push_back( xIdx_[inst_.machine(j,h-1)][j][t] );
                coef.push_back( t + inst_.time(j,inst_.machine(j,h-1)) );
                idx.push_back( xIdx_[inst_.machine(j,h)][j][t] );
                coef.push_back( -1 * t );
            }

            lp_add_row( mip, idx, coef, "ord("+to_string(h+1)+","+to_string(j+1)+")", 'L', 0 );
            
        }
    }

    cout << "ord restriction added" << endl;
    lp_write_lp( mip, inst_.instanceName().c_str() );
    lp_write_mps( mip, inst_.instanceName().c_str() );
    if (inst_.execute()){
        lp_optimize( mip );
        lp_write_sol(mip, "jssp_Kondili.sol");
    }
}


Kondili::~Kondili()
{
    lp_free( &mip );
}

