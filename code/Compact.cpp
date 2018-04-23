#include "Compact.hpp"

#include <vector>
#include <string>
#include <cfloat>

using namespace std;

Compact::Compact( const Instance &_inst ) : inst_(_inst) { // já inicializa a variável privada _inst com o valor passado por referencia
    mip = lp_create();

    // variáveis de decisão
    xIdx_ = vector<vector<int>>(inst_.n(),vector<int>(inst_.m()));
    yIdx_ = vector<vector<vector<int>>>(inst_.n(),vector<vector<int>>(inst_.n(),vector<int>(inst_.m(),-1))); // inicia todos os valores com -1

    vector< string > names; // nome das variáveis
    vector< double > lb; // lower bound
    vector< double > ub; // upper bound
    vector< double > obj; // se é objetivo?
    vector< char > integer; // variável inteira?

    // criação das variáveis x
    for (int j = 0; j < inst_.n(); j++){
        for (int a = 0; a < inst_.m(); a++){
            xIdx_[j][a] = names.size(); // número do índice da variável (vai de 1 até n*m)
            names.push_back("x("+to_string(j+1)+","+to_string(a+1)+")"); // nome dessa variável
            lb.push_back(0.0);
            ub.push_back(DBL_MAX);
            obj.push_back(0.0);
            integer.push_back(1);
        }
    }
    // criação das variáveis y
    for (int j = 0; j < inst_.n(); j++){
        for (int i = 0; i < inst_.n(); i++){
            for (int a = 0; a < inst_.m(); a++){
                yIdx_[j][i][a] = names.size();
                names.push_back("y("+to_string(j+1)+","+to_string(i+1)+","+to_string(a+1)+")");
                lb.push_back(0.0);
                ub.push_back(1.0); // variável binária
                obj.push_back(0.0);
                integer.push_back(1);
            }
        }
    }

    // c var
    cIdx_ = names.size();
    names.push_back("Z");
    lb.push_back( 0.0 );
    ub.push_back( DBL_MAX );
    obj.push_back( 1.0 );
    integer.push_back( 1 );

    // adiciona colunas ao solver
    lp_add_cols( mip, obj, lb, ub, integer, names );

    // restrição ord (cada máquina só pode processar após a anterior ter sido processada. x[j,sigma[j,a]] - x[j,sigma[j,a-1]] >= p[j][a-1]
    for (int j = 0; j < inst_.n(); j++){
        for (int a = 1; a < inst_.m(); a++){
            vector< int > idx;
            vector< double > coef;

            idx.push_back( xIdx_[j][inst_.machine(j,a)] );
            coef.push_back( 1.0 );
            idx.push_back( xIdx_[j][inst_.machine(j,a-1)] );
            coef.push_back( -1.0 );

            // adiciona restrição.
            lp_add_row( mip, idx, coef, "ord("+to_string(j+1)+","+to_string(a+1)+")", 'G', inst_.time( j, inst_.machine(j,a-1)) );
        }
    }

    double K; // variável INF
    for (int j = 0; j < inst_.n(); j++){
        for (int a = 0; a < inst_.m(); a++){
            K += inst_.time(j,a);
        }
    }

    // restrição phi e psi

    for (int i = 0; i < inst_.n(); i++){
        for (int j = 0; j < inst_.n(); j++){
            if (i == j) continue;
            for (int a = 0; a < inst_.m(); a++){
                vector< int > idx_phy;
                vector< double > coef_phy;

                idx_phy.push_back(xIdx_[i][a]);
                coef_phy.push_back(1.0);
                idx_phy.push_back(xIdx_[j][a]);
                coef_phy.push_back(-1.0);
                idx_phy.push_back(yIdx_[i][j][a]);
                coef_phy.push_back(K);

                lp_add_row( mip, idx_phy, coef_phy, "phi("+to_string(i+1)+","+to_string(j+1)+","+to_string(a+1)+")", 'G', inst_.time( j, a) );

                vector< int > idx_psi;
                vector< double > coef_psi;

                idx_psi.push_back(xIdx_[j][a]);
                coef_psi.push_back(1.0);
                idx_psi.push_back(xIdx_[i][a]);
                coef_psi.push_back(-1.0);
                idx_psi.push_back(yIdx_[i][j][a]);
                coef_psi.push_back(-1.0*K);
                double c = inst_.time(i,a) - K;

                lp_add_row( mip, idx_psi, coef_psi, "psi("+to_string(i+1)+","+to_string(j+1)+","+to_string(a+1)+")", 'G', c );
            }
        }
    }

    // restrições fin
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        idx.push_back( cIdx_ );
        coef.push_back( 1.0 );

        idx.push_back( xIdx_[j][inst_.machine(j, inst_.m()-1)] );
        coef.push_back( -1.0 );

        lp_add_row( mip, idx, coef, "fin("+to_string(j+1)+")", 'G', inst_.time( j, inst_.machine(j, inst_.m()-1)) );
    }

    lp_optimize( mip );
    lp_write_lp( mip, "jssp_compact" );
    lp_write_sol(mip, "jssp_compact.sol");
}

// Compact::Compact( const Instance &_inst ) :
//     inst_(_inst)
// {
//     mip = lp_create();

//     vector< string > names;
//     vector< double > lb;
//     vector< double > ub;
//     vector< double > obj;
//     vector< char > integer;

//     xIdx_ = vector< vector< int > >( inst_.n(), vector<int>( inst_.m(), -1 ) );

//     yIdx_ = vector< vector< vector<int> > >( inst_.n(), vector< vector< int > >( inst_.n(),
//                 vector<int>( inst_.m(), -1 ) ) );

//     // creating x vars
//     for ( int j=0 ; (j<inst_.n()) ; ++j )
//     {
//         for ( int i=0 ; (i<inst_.m()) ; ++i )
//         {
//             xIdx_[j][i] = names.size();
//             names.push_back( "x("+to_string(j+1)+","+to_string(i+1)+")" );
//             lb.push_back( 0.0 );
//             ub.push_back( inst_.lst( j, i ) );
//             obj.push_back( 0.0 );
//             integer.push_back( 1 );
//         }
//     }

//     // y vars
//     for ( int j1=0 ; (j1<inst_.n()) ; ++j1 )
//     {
//         for ( int j2=j1+1 ; (j2<inst_.n()) ; ++j2 )
//         {
//             for ( int i=0 ; (i<inst_.m()) ; ++i )
//             {
//                 yIdx_[j1][j2][i] = names.size();
//                 names.push_back( "y("+to_string(j1+1)+","+to_string(j2+1)+","+to_string(i+1)+")" );
//                 lb.push_back( 0.0 );
//                 ub.push_back( 1.0 );
//                 obj.push_back( 0.0 );
//                 integer.push_back( 1 );
//              }
//         }
//     }

//     // c var
//     cIdx_ = names.size();
//     names.push_back("C");
//     lb.push_back( 0.0 );
//     ub.push_back( DBL_MAX );
//     obj.push_back( 1.0 );
//     integer.push_back( 1 );

//     lp_add_cols( mip, obj, lb, ub, integer, names );

//     // constraint for job on each one of its machines
//     for ( int j=0 ; (j<inst_.n()) ; ++j )
//     {
//         for ( int i=1 ; (i<inst_.m()) ; ++i )
//         {
//             vector< int > idx;
//             vector< double > coef;

//             idx.push_back( xIdx_[j][inst_.machine(j,i)] );
//             coef.push_back( 1.0 );
//             idx.push_back( xIdx_[j][inst_.machine(j,i-1)] );
//             coef.push_back( -1.0 );

//             lp_add_row( mip, idx, coef, "prec("+to_string(j+1)+","+to_string(i+1)+")", 'G', inst_.time( j, inst_.machine(j,i-1)) );
//         }
//     }

//     // linking c and x
//     for ( int j=0 ; (j<inst_.n()) ; ++j )
//     {
//         vector< int > idx;
//         vector< double > coef;

//         idx.push_back( cIdx_ );
//         coef.push_back( 1.0 );

//         idx.push_back( xIdx_[j][inst_.machine(j, inst_.m()-1)] );
//         coef.push_back( -1.0 );

//         lp_add_row( mip, idx, coef, "lnkCX("+to_string(j+1)+")", 'G', inst_.time( j, inst_.machine(j, inst_.m()-1)) );
//     }

//     for ( int j1=0 ; (j1<inst_.n()) ; ++j1 )
//     {
//         for ( int j2=0 ; (j2<inst_.n()) ; ++j2 )
//         {
//             if (j1==j2)
//                 continue;
            
//             for ( int i=0 ; (i<inst_.m()) ; ++i )
//             {
//                 vector< int > idx; vector< double > coef;
//                 idx.push_back( xIdx_[j1][i] );
//                 coef.push_back( 1.0 );
//                 idx.push_back( xIdx_[j2][i] );
//                 coef.push_back( -1.0 );
    
//                 double rhs = inst_.time( j2, i );

//                 if (j1<j2)
//                 {
//                     rhs -= 99999;
//                     idx.push_back( yIdx_[j1][j2][i] );
//                     coef.push_back( -99999 );
//                 }
//                 else
//                 {
//                     idx.push_back( yIdx_[j2][j1][i] );
//                     coef.push_back( 99999 );
//                 }

//                 lp_add_row( mip, idx, coef, "lnkXY("+to_string(j1+1)+","+to_string(j2+1)+","+to_string(i+1)+")", 'G', rhs );
//             }
//         }
//     }
        
//     lp_optimize( mip );
//     lp_write_lp( mip, "jssp" );
// }

Compact::~Compact()
{
    lp_free( &mip );
}

