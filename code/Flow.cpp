#include "Flow.hpp"

#include <vector>
#include <string>
#include <cfloat>
#include <iostream>
#include <fstream>

using namespace std;

Flow::Flow( const Instance &_inst ) :
    inst_(_inst) {
    mip = lp_create();

    vector< string > names;
    vector< double > lb;
    vector< double > ub;
    vector< double > obj;
    vector< char > integer;

    xIdx_ = vector<vector<vector<vector<vector<int>>>>> (inst_.n(), vector<vector<vector<vector<int>>>>(inst_.m()+1, vector<vector<vector<int>>>(inst_.maxTime(),vector<vector<int>>(inst_.m()+2, vector<int>(inst_.maxTime()))))) ; 
    

    // creating x vars
    for ( int j=0 ; (j<inst_.n()) ; ++j ) {
        for ( int m0=0 ; (m0<=inst_.m()) ; ++m0 ) {
            for ( int t0 = 0; t0 < inst_.maxTime()-1; t0++  ) {
                for (int mf = 1; mf <= inst_.m()+1; mf++){
                    for (int tf=t0+1; tf < inst_.maxTime(); tf++){
                        xIdx_[j][m0][t0][mf][tf] = names.size();
                        if (m0 == mf && tf > t0+1) continue;
                        if ( (m0 == 0 && mf == inst_.m()+1) )  continue; // variáveis que não entram
                        if (m0 == 0){
                            names.push_back( "x("+to_string(j+1)+",i,"+to_string(t0+1)+","+to_string(mf)+","+to_string(tf+1)+")" );
                        } else
                        if (mf == inst_.m()+1){
                            names.push_back( "x("+to_string(j+1)+","+to_string(m0)+","+to_string(t0+1)+",f,"+to_string(tf+1)+")" );
                        } else {
                            names.push_back( "x("+to_string(j+1)+","+to_string(m0)+","+to_string(t0+1)+","+to_string(mf)+","+to_string(tf+1)+")" );
                        }
                        lb.push_back( 0.0 );
                        ub.push_back( 1 );
                        obj.push_back( 0 );
                        integer.push_back( 1 );
                    }
                }
            }
        }
    }
    ofstream f;
    f.open ("variables.txt");
    for (string name : names){
        f << name << endl;
    }
    f.close();
    
    cout << "variáveis criadas" << endl;

    // c var
    cIdx_ = names.size();
    names.push_back("C");
    lb.push_back( 0.0 );
    ub.push_back( DBL_MAX );
    obj.push_back( 1.0 );
    integer.push_back( 1 );

    //lp_add_cols( mip, obj, lb, ub, integer, names );

    // set of entering flows
    vector<vector<vector<int>>> enter_flow = vector<vector<vector<int>>>(inst_.m()+2,vector<vector<int>>(inst_.maxTime()));
    // set of exiting flows
    vector<vector<vector<int>>> exit_flow = vector<vector<vector<int>>>(inst_.m()+1,vector<vector<int>>(inst_.maxTime()));

    for ( int j=0 ; (j<inst_.n()) ; ++j ) {
        for ( int m0=0 ; (m0<=inst_.m()) ; ++m0 ) {
            for ( int t0 = 0; t0 < inst_.maxTime()-1; t0++  ) {
                for (int mf = 1; mf <= inst_.m()+1; mf++){
                    for (int tf=t0+1; tf < inst_.maxTime(); tf++){
                        enter_flow[mf][tf].push_back(xIdx_[j][m0][t0][mf][tf]);
                        exit_flow[m0][t0].push_back(xIdx_[j][m0][t0][mf][tf]);
                    }
                }
            }
        }
    }

    cout << "fluxos criados" << endl;

    f.open ("exit_flows.txt");
    for (int mf = 1; mf < inst_.m()+1; mf++){
        for (int tf=1; tf < inst_.maxTime(); tf++){
            f << "machine " << mf << " time " << tf << endl;
            for (int var : exit_flow[mf][tf]){
                f << names[var] << endl;
            }
        }
    }
    f.close();
    cout << "exit_flows criado" << endl;

    f.open ("enter_flows.txt");
    for (int m0 = 0; m0 <= inst_.m(); m0++){
        for (int t0=0; t0 < inst_.maxTime()-1; t0++){
            f << "machine " << m0 << " time " << t0 << endl;
            for (int var : enter_flow[m0][t0]){
                f << names[var] << endl;
            }
        }
    }
    f.close();
    cout << "enter_flows criado" << endl;
    createCompleteGraphDot();



    cout << "variáveis criadas" << endl;

    // constraint for job on each one of its machines
    // for ( int j=0 ; (j<inst_.n()) ; ++j )
    // {
    //     for ( int i=1 ; (i<inst_.m()) ; ++i )
    //     {
    //         vector< int > idx;
    //         vector< double > coef;

    //         idx.push_back( xIdx_[j][inst_.machine(j,i)] );
    //         coef.push_back( 1.0 );
    //         idx.push_back( xIdx_[j][inst_.machine(j,i-1)] );
    //         coef.push_back( -1.0 );

    //         lp_add_row( mip, idx, coef, "prec("+to_string(j+1)+","+to_string(i+1)+")", 'G', inst_.time( j, inst_.machine(j,i-1)) );
    //     }
    // }

    // // linking c and x
    // for ( int j=0 ; (j<inst_.n()) ; ++j )
    // {
    //     vector< int > idx;
    //     vector< double > coef;

    //     idx.push_back( cIdx_ );
    //     coef.push_back( 1.0 );

    //     idx.push_back( xIdx_[j][inst_.machine(j, inst_.m()-1)] );
    //     coef.push_back( -1.0 );

    //     lp_add_row( mip, idx, coef, "lnkCX("+to_string(j+1)+")", 'G', inst_.time( j, inst_.machine(j, inst_.m()-1)) );
    // }

    // for ( int j1=0 ; (j1<inst_.n()) ; ++j1 )
    // {
    //     for ( int j2=0 ; (j2<inst_.n()) ; ++j2 )
    //     {
    //         if (j1==j2)
    //             continue;
            
    //         for ( int i=0 ; (i<inst_.m()) ; ++i )
    //         {
    //             vector< int > idx; vector< double > coef;
    //             idx.push_back( xIdx_[j1][i] );
    //             coef.push_back( 1.0 );
    //             idx.push_back( xIdx_[j2][i] );
    //             coef.push_back( -1.0 );
    
    //             double rhs = inst_.time( j2, i );

    //             if (j1<j2)
    //             {
    //                 rhs -= 99999;
    //                 idx.push_back( yIdx_[j1][j2][i] );
    //                 coef.push_back( -99999 );
    //             }
    //             else
    //             {
    //                 idx.push_back( yIdx_[j2][j1][i] );
    //                 coef.push_back( 99999 );
    //             }

    //             lp_add_row( mip, idx, coef, "lnkXY("+to_string(j1+1)+","+to_string(j2+1)+","+to_string(i+1)+")", 'G', rhs );
    //         }
    //     }
    // }
        
    // lp_optimize( mip );
    // lp_write_lp( mip, "jssp" );
}

Flow::~Flow()
{
    lp_free( &mip );
}

// void Flow::create_x11(){
//     x11_color.push_back("black");
//     x11_color.push_back("magenta");
//     x11_color.push_back("grey");
//     x11_color.push_back("gold");
//     x11_color.push_back("blue");
//     x11_color.push_back("green");
//     x11_color.push_back("yellow");
//     x11_color.push_back("red");
//     x11_color.push_back("cyan");
//     x11_color.push_back("crimson");
//     x11_color.push_back("chocolate");
//     x11_color.push_back("brown");
//     x11_color.push_back("orange");
// }

void Flow::createCompleteGraphDot(){
    ofstream f;
    f.open ("complete.dot");

    f << "digraph complete {" << endl;
    f << "ratio = \"auto\" ;" << endl;
    f << "rankdir=LR;" << endl << endl;
    for (int m = 0; m <= inst_.m()+1; m++){
        string mach = to_string(m);
        if (mach == "0") mach = "i";
        else if (mach == to_string(inst_.m()+1)) mach = "f";
        for (int t = 0; t < inst_.maxTime(); t++){
            f << "\"" << mach << "," << t+1 << "\" [shape=box, regular=1,style=filled,fillcolor=white];" << endl;
        }
    }
    string mach0,machf;
    for ( int j=0 ; (j<inst_.n()) ; ++j ) {
        for ( int m0=0 ; (m0<=inst_.m()) ; ++m0 ) {
            mach0 = to_string(m0);
            if (mach0 == "0") mach0 = "i";
            for ( int t0 = 0; t0 < inst_.maxTime()-1; t0++  ) {
                for (int mf = 1; mf <= inst_.m()+1; mf++){
                    machf = to_string(mf);
                    if (machf == to_string(inst_.m()+1)) machf = "f";
                    for (int tf=t0+1; tf < inst_.maxTime(); tf++){
                        if (m0 == mf && tf > t0+1) continue;
                        f << "\"" << m0 << "," << t0+1 << "\" -> " << "\"" << mf << "," << tf << "\"" << "[ label=\"x(" << j+1 << "," << mach0 << "," << t0+1 << "," << machf << "," << tf+1 << ")\" ];"<< endl ; // \" color=\"" << x11_color[j] <<
                    }
                }
            }
        }
    }

    for (int m = 0; m <= inst_.m()+1; m++){
        string mach = to_string(m);
        if (mach == "0") mach = "i";
        else if (mach == to_string(inst_.m()+1)) mach = "f";
        f << "{rank=same;";
        for (int t = 0; t < inst_.maxTime(); t++){
            f << " \"" << mach << "," << t+1 << "\"";
        }
        f << "}" << endl;
    }
    f << "}" << endl;
    f.close();
}