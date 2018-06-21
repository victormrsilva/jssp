#include "Fernando.hpp"

#include <vector>
#include <string>
#include <cfloat>
#include <iostream>
#include <fstream>
using namespace std;

Fernando::Fernando( const Instance &_inst ) : inst_(_inst) { // já inicializa a variável privada _inst com o valor passado por referencia
    mip = lp_create();

    // variáveis de decisão
    xIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    eIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    fIdx_ = vector<vector<int>>(inst_.m(),vector<int>(inst_.maxTime()+1));
    

    vector< string > names; // nome das variáveis
    vector< double > lb; // lower bound
    vector< double > ub; // upper bound
    vector< double > obj; // se é objetivo?
    vector< char > integer; // variável inteira?

    // criação das variáveis x
    for (int i = 0; i < inst_.m(); i++){
        for (int t = 0; t < inst_.maxTime(); t++){
            fIdx_[i][t] = names.size(); 
            names.push_back("f("+to_string(i+1)+","+to_string(t+1)+")"); // nome dessa variável
            lb.push_back(0.0);
            ub.push_back(1.0);
            obj.push_back(0.0);
            integer.push_back(1);
    
            for (int j = 0; j < inst_.n(); j++){

                if (t >= inst_.est(j,i) && t <= inst_.lst(j,i)){
                    xIdx_[i][j][t] = names.size(); 
                    names.push_back("x("+to_string(i+1)+","+to_string(j+1)+","+to_string(t+1)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(1);

                
                    eIdx_[i][j][t] = names.size(); 
                    names.push_back("e("+to_string(i+1)+","+to_string(j+1)+","+to_string(t+1)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(1);
                }
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

    ofstream f;
    f.open ("variables.txt");
    for (string name : names){
        f << name << endl;
    }
    f.close();

    // adiciona colunas ao solver
    lp_add_cols( mip, obj, lb, ub, integer, names );
    cout << "Number of variables: " << names.size() << endl;
    //restriction 28
    for (int i = 0; i < inst_.m(); i++){
        vector< int > idx;
        vector< double > coef;

        for (int j = 0; j < inst_.n(); j++){
            if (inst_.machine(j,0) == i) { // se for a primeira máquina
                idx.push_back( xIdx_[i][j][0] );
                coef.push_back( -1.0 );
            }
            // adiciona restrição.
        }
        idx.push_back( fIdx_[i][0] );
        coef.push_back( -1.0 );

        lp_add_row( mip, idx, coef, "inicio(m"+to_string(i+1)+",t"+to_string(1)+")", 'E', -1 );
    }
    cout << "restriction 28 added" << endl;
    // restriction 29
    for (int i = 0; i < inst_.m(); i++){
        for (int t = 1; t < inst_.maxTime(); t++){
            vector< int > idx;
            vector< double > coef;

            for (int j = 0; j < inst_.n(); j++){
                int tp = t - inst_.time(j,i);
                //cout << "m = " << i << " t = " << t << " j = " << j << " tp = " << tp << " p = " << inst_.time(j,i) << " est " << inst_.est(j,i) << " lst " << inst_.lst(j,i) << endl;
                if (tp >=  inst_.est(j,i)  && tp <= inst_.lst(j,i)){
                    idx.push_back( xIdx_[i][j][tp] );
                    coef.push_back( 1.0 );
                }
                if (t >= inst_.est(j,i) && t <= inst_.lst(j,i)){
                    idx.push_back( xIdx_[i][j][t] );
                    coef.push_back( -1.0 );
                }
                // adiciona restrição.
            }
            idx.push_back( fIdx_[i][t-1] );
            coef.push_back( 1.0 );
            idx.push_back( fIdx_[i][t] );
            coef.push_back( -1.0 );

            lp_add_row( mip, idx, coef, "c29("+to_string(i+1)+","+to_string(t+1)+")", 'E', 0 );
        }

    }
    cout << "restriction 29 added" << endl;

    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        int h = inst_.machine(j,0); // first machine

        idx.push_back( xIdx_[h][j][0] );
        coef.push_back( -1.0 );
        // adiciona restrição.
        idx.push_back( eIdx_[h][j][0] );
        coef.push_back( -1.0 );

        lp_add_row( mip, idx, coef, "espera_ini(m"+to_string(h)+",j"+to_string(j+1)+",t"+to_string(1)+")", 'E', -1 );

    }
    cout << "restriction 30 added" << endl;

    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            int h = inst_.machine(j,i); // first machine

            for (int t = inst_.est(j,h); t < inst_.lst(j,h); t++){
                if (t == 0) continue;
                vector< int > idx;
                vector< double > coef;
                if (i != 0){
                    int h_anterior = inst_.machine(j,i-1); // first machine

                    int tp = t - inst_.time(j,h_anterior);

                    if (tp >=  inst_.est(j,h_anterior)  && tp <= inst_.lst(j,h_anterior)){
                        idx.push_back( xIdx_[h_anterior][j][tp] );
                        coef.push_back( 1.0 );
                    }
                }
                if (t > inst_.est(j,h)){
                    idx.push_back( eIdx_[h][j][t-1] );
                    coef.push_back( 1.0 );
                }

                idx.push_back( xIdx_[h][j][t] );
                coef.push_back( -1.0 );
                // adiciona restrição.
                idx.push_back( eIdx_[h][j][t] );
                coef.push_back( -1.0 );

                lp_add_row( mip, idx, coef, "execute_wait("+to_string(h+1)+","+to_string(j+1)+","+to_string(t+1)+")", 'E', 0 );

            }
        }
    }
    cout << "restriction 31 added" << endl;

    // for (int i = 1; i < inst_.m(); i++){
    //     for (int j = 0; j < inst_.n(); j++){
    //         int h = inst_.machine(j,i); // first machine
    //         int h_anterior = inst_.machine(j,i-1); // first machine
    //         for (int t = inst_.est(j,h); t < inst_.lst(j,h); t++){
    //             vector< int > idx;
    //             vector< double > coef;


    //             int t_h_anterior = t - inst_.time(j,h_anterior);
    //             if (t_h_anterior >= inst_.est(j,h_anterior) && t_h_anterior  <= inst_.lst(j,h_anterior) ){
    //                 idx.push_back( xIdx_[h_anterior][j][t - inst_.time(j,h_anterior)] );
    //                 coef.push_back( 1.0 );
    //             }

    //             if (t > 0){
    //                 idx.push_back( eIdx_[h][j][t-1] );
    //                 coef.push_back( 1.0 );
    //             }

    //             idx.push_back( xIdx_[h][j][t] );
    //             coef.push_back( -1.0 );
    //             // adiciona restrição.
    //             idx.push_back( eIdx_[h][j][t] );
    //             coef.push_back( -1.0 );

    //             lp_add_row( mip, idx, coef, "c32("+to_string(i+1)+","+to_string(j+1)+","+to_string(t+1)+")", 'E', 0 );

    //         }
    //     }
    // }
    // cout << "restriction 32 added" << endl;

    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            vector< int > idx;
            vector< double > coef;
            int h = inst_.machine(j,i);
            int t = inst_.lst(j,h);
            idx.push_back( xIdx_[h][j][t] );
            coef.push_back( -1.0 );
            idx.push_back( eIdx_[h][j][t-1] );
            coef.push_back( 1.0 );
            if (i > 0) {
                int h_anterior = inst_.machine(j,i-1);
                int tp = t - inst_.time(j,h_anterior);
                idx.push_back( xIdx_[h_anterior][j][tp] );
                coef.push_back( 1.0 );
            }
            lp_add_row( mip, idx, coef, "final_execute("+to_string(h+1)+","+to_string(j+1)+")", 'E', 0.0 );
        }
    }
    cout << "final_execute constraints created" << endl;

    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        idx.push_back( cIdx_ );
        coef.push_back( 1.0 );
        int h = inst_.machine(j,inst_.m()-1);
        for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
            
            idx.push_back(xIdx_[h][j][t]);
            coef.push_back(-1*(t+inst_.time(j,h)));
        }

        lp_add_row( mip, idx, coef, "makespan("+to_string(j+1)+")", 'G', 0 );
    }
    cout << "end constraints created" << endl;

    lp_write_lp( mip, inst_.instanceName().c_str() );
    lp_write_mps( mip, inst_.instanceName().c_str() );
    if (inst_.execute()){
        lp_optimize( mip );
        lp_write_sol(mip, "jssp_Fernando.sol");
    }
}


Fernando::~Fernando()
{
    lp_free( &mip );
}

