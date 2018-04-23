#include "Flow.hpp"

#include <vector>
#include <string>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

int Flow::getXidx(int j, int m0, int t0, int mf, int tf) const{
    auto it = xIdx_[j][m0].find(t0);
    if (it == xIdx_[j][m0].end()) return -1;
    auto it2 = it->second.find(mf);
    if (it2 == it->second.end()) return -1;
    auto it3 = it2->second.find(tf);
    if (it3 == it2->second.end()) return -1;
    return it3->second; // return value of map tf (number of variable)
}

Flow::Flow( const Instance &_inst ) :
    inst_(_inst),
    mip(lp_create()),
    xIdx_(vector< vector< map< int, map< int, map< int, int >>>>> (inst_.n(), vector< map< int, map< int, map< int, int >>>>(inst_.m()+1))),
    enter_flow(vector<vector<vector<vector<int>>>>(inst_.n(),vector<vector<vector<int>>>(inst_.m()+2+inst_.n(),vector<vector<int>>(inst_.maxTime())))),
    exit_flow(vector<vector<vector<vector<int>>>>(inst_.n(),vector<vector<vector<int>>>(inst_.m()+1,vector<vector<int>>(inst_.maxTime())))),
    process(vector<vector<vector<vector<int>>>>(inst_.n(),(vector<vector<vector<int>>>(inst_.m()+1,vector<vector<int>>(inst_.maxTime()))))) {
    
    vector< string > names;
    vector< double > lb;
    vector< double > ub;
    vector< double > obj;
    vector< char > integer;

    // creating x vars
    for ( int j=0 ; (j<inst_.n()) ; ++j ) {
        for ( int m=-1 ; (m < inst_.m()) ; ++m ) {
            if (m == -1){ // máquina inicial
                
                xIdx_[j][m+1][0][inst_.machine(j,0)+1][0] = names.size();
                exit_flow[j][0][j].push_back(names.size());
                enter_flow[j][inst_.machine(j,0)+1][0].push_back(names.size());
                names.push_back( "x("+to_string(j+1)+",i,0,"+to_string(inst_.machine(j,0)+1)+",0)" );
                lb.push_back( 0.0 );
                ub.push_back( 1 );
                obj.push_back( 0 );
                integer.push_back( 1 );
                continue;
            }
            int m0 = inst_.machine(j,m);
            int mf = (m == inst_.m()-1 ? inst_.m() : inst_.machine(j,m+1));
            int dur = inst_.time(j,m0); // duration time for machine m0 in job j
            for (int t = inst_.est(j,m0); t < inst_.lst(j,m0); t++){
                if (mf < inst_.m()){
                    // arc for another machine
                    xIdx_[j][m0+1][t][mf+1][t+dur] = names.size();
                    //cout << j << " " << m0+1 << " " << t << " " << mf+1 << " " << t+dur << endl;
                    enter_flow[j][mf+1][t+dur].push_back(names.size());
                    exit_flow[j][m0+1][t].push_back(names.size());
                    for (int tp = t; tp < t+dur; tp++){
                        process[j][m0+1][tp].push_back(names.size());
                    }
                    names.push_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(mf+1)+","+to_string(t+dur)+")" );
                    lb.push_back( 0.0 );
                    ub.push_back( 1 );
                    obj.push_back( 0 );
                    integer.push_back( 1 );
                    // arc for same machine (waiting) 
                    // can only be made in the last moment possible
                    if (t == inst_.lst(j,m0)-1) continue;
                    // else
                    xIdx_[j][m0+1][t][m0+1][t+1] = names.size();
                    enter_flow[j][m0+1][t+1].push_back(names.size());
                    exit_flow[j][m0+1][t].push_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.push_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(m0+1)+","+to_string(t+1)+")" );
                    lb.push_back( 0.0 );
                    ub.push_back( 1 );
                    obj.push_back( 0 );
                    integer.push_back( 1 );
                } else { // conclusion machine f
                    xIdx_[j][m0+1][t][mf+1][t+dur] = names.size();
                    enter_flow[j][mf+1+j][t+dur].push_back(names.size());
                    exit_flow[j][m0+1][t].push_back(names.size());
                    for (int tp = t; tp < t+dur; tp++){
                        process[j][m0+1][tp].push_back(names.size());
                    }
                    names.push_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+",f,"+to_string(t+dur)+")" );
                    lb.push_back( 0.0 );
                    ub.push_back( 1 );
                    obj.push_back( 0 );
                    integer.push_back( 1 );
                    
                    xIdx_[j][m0+1][t][m0+1][t+1] = names.size();
                    enter_flow[j][m0+1][t+1].push_back(names.size());
                    exit_flow[j][m0+1][t].push_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.push_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(m0+1)+","+to_string(t+1)+")" );
                    lb.push_back( 0.0 );
                    ub.push_back( 1 );
                    obj.push_back( 0 );
                    integer.push_back( 1 );
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

    lp_add_cols( mip, obj, lb, ub, integer, names );

    f.open ("exit_flows.txt");
    for (int j = 0; j < inst_.n(); j++){
        for (int tf=0; tf < inst_.maxTime(); tf++){
            for (int mf = 0; mf < inst_.m()+1; mf++){
            
                f << "machine " << mf << " time " << tf << endl;
                for (int var : exit_flow[j][mf][tf]){
                    f << names[var] << endl;
                }
            }
        }
    }
    f.close();
    cout << "exit_flows criado" << endl;

    f.open ("enter_flows.txt");
    for (int j = 0; j < inst_.n(); j++){
        for (int t0=0; t0 < inst_.maxTime()-1; t0++){
            for (int m0 = 0; m0 <= inst_.m()+1; m0++){
            
                f << "machine " << m0 << " time " << t0 << endl;
                for (int var : enter_flow[j][m0][t0]){
                    f << names[var] << endl;
                }
            }
        }
    }
    f.close();
    cout << "enter_flows criado" << endl;

    f.open ("processing_machines.txt");
    for (int t0=0; t0 < inst_.maxTime()-1; t0++){
        for (int m0 = 0; m0 <= inst_.m(); m0++){
        
            f << "machine " << m0 << " time " << t0 << endl;
            for (int j = 0; j < inst_.n(); j++){
                for (int var : process[j][m0][t0]){
                    f << names[var] << endl;
                }
            }
        }
    }
    f.close();
    cout << "processing_machines created" << endl;
    createCompleteGraphDot();
    cout << "arquivo dot criado" << endl;

    // initial flow constraint
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;
        for (int var : exit_flow[j][0][j]){
            idx.push_back( var );
            coef.push_back( -1.0 );
        }
        lp_add_row( mip, idx, coef, "init_flow("+to_string(j)+")", 'E', -1.0 );
    }
    cout << "initial flow constraints ok" << endl;

    // final flow constraint
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;
        for (int t = 1; t < inst_.maxTime(); t++){
            for (int var : enter_flow[j][inst_.m()+1+j][t]){
                idx.push_back( var );
                coef.push_back( 1.0 );
            }
        }

        if (idx.size() != 0){
            lp_add_row( mip, idx, coef, "final_flow("+to_string(j)+")", 'E', 1.0 );
        }
    }
    cout << "final flow constraints ok" << endl;
    // flow constraints
    for (int j = 0; j < inst_.n(); j++){
        for ( int t=0 ; (t<inst_.maxTime()) ; ++t ) {
            for ( int m=1 ; (m<=inst_.m()) ; ++m ) {
                vector< int > idx;
                vector< double > coef;
                for (int var : enter_flow[j][m][t]){
                    idx.push_back( var );
                    coef.push_back( 1.0 );
                }
                for (int var : exit_flow[j][m][t]){
                    idx.push_back( var );
                    coef.push_back( -1.0 );
                }

                if (idx.size() != 0){
                    lp_add_row( mip, idx, coef, "flow("+to_string(j+1)+","+to_string(m)+","+to_string(t)+")", 'E', 0.0 );
                }
            }
        }
    }
    cout << "flow constraints created" << endl;

    // processing restriction
    for ( int t=0 ; (t<inst_.maxTime()) ; ++t ) {

        for ( int m=0 ; (m<=inst_.m()) ; ++m ) {
            vector< int > idx;
            vector< double > coef;
            for (int j = 0; j < inst_.n(); j++){
                for( int var : process[j][m][t]){
                    idx.push_back( var );
                    coef.push_back( 1.0 );

                }
            }
            if (idx.size() != 0){
                lp_add_row( mip, idx, coef, "processing("+to_string(m) + ","+to_string(t) + ")", 'L', 1.0 );        
            }
        }
        
    }
    cout << "processing constraints created" << endl;

    // restrições fin
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        idx.push_back( cIdx_ );
        coef.push_back( 1.0 );

        for (int t = 1; t < inst_.maxTime(); t++){
            
            for (int var : enter_flow[j][inst_.m()+1+j][t]){
                idx.push_back( var );
                coef.push_back( -t );
            }
        }

        lp_add_row( mip, idx, coef, "fin("+to_string(j+1)+")", 'G', 0 );
    }
    cout << "end constraints created" << endl;

// restrições opt
    vector< int > idx;
    vector< double > coef;

    idx.push_back( xIdx_[0][1][6][2][6+inst_.time(0,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[0][2][25][4][25+inst_.time(0,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[0][3][5][1][5+inst_.time(0,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[0][4][31][6][31+inst_.time(0,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[0][5][42][7][42+inst_.time(0,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[0][6][38][5][38+inst_.time(0,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(1,6)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][1][38][4][38+inst_.time(1,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][2][0][3][0+inst_.time(1,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][3][8][5][8+inst_.time(1,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][4][48][7][48+inst_.time(1,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][5][15][6][15+inst_.time(1,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[1][6][28][1][28+inst_.time(1,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(2,6)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][1][19][2][19+inst_.time(2,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][2][31][5][31+inst_.time(2,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][3][0][4][0+inst_.time(2,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][4][7][6][7+inst_.time(2,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][5][48][7][48+inst_.time(2,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[2][6][11][1][11+inst_.time(2,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(3,6)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][1][13][3][13+inst_.time(3,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][2][8][1][8+inst_.time(3,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][3][22][4][22+inst_.time(3,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][4][27][5][27+inst_.time(3,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][5][30][6][30+inst_.time(3,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[3][6][46][7][46+inst_.time(3,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(4,6)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][1][48][4][48+inst_.time(4,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][2][22][5][22+inst_.time(4,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][3][13][2][13+inst_.time(4,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][4][52][7][52+inst_.time(4,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][5][25][6][25+inst_.time(4,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[4][6][41][1][41+inst_.time(4,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(5,6)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][1][28][5][28+inst_.time(5,0)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,1)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][2][13][4][13+inst_.time(5,1)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,2)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][3][54][7][54+inst_.time(5,2)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,3)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][4][16][6][16+inst_.time(5,3)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,4)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][5][38][3][38+inst_.time(5,4)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,5)", 'E', 1 );
    idx.clear();
    coef.clear();

    idx.push_back( xIdx_[5][6][19][1][19+inst_.time(5,5)] );
    coef.push_back( 1.0 );
    lp_add_row( mip, idx, coef, "opt(6,6)", 'E', 1 );
    idx.clear();
    coef.clear();


    cout << "optmium constraints ok" << endl;
        
     lp_optimize( mip );
     lp_write_lp( mip, "jssp" );
     lp_write_sol(mip, "jssp.sol");
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
    //f << "ratio = \"auto\" ;" << endl;
    f << "rankdir=LR;" << endl << endl;
    for (int m = 0; m <= inst_.m()+1; m++){
        if (m == 0){
            f << "\"i,1\" [shape=box];" << endl;
            continue;
        }
        string mach;
        if (mach == to_string(inst_.m()+1)) mach = "f";
        else mach = to_string(m);
        for (int t = 0; t < inst_.maxTime(); t++){
            f << "\"" << mach << "," << t+1 << "\" [shape=box];" << endl;
        }
    }

    for ( int j=0 ; (j<inst_.n()) ; ++j ) {
        for ( int m=-1 ; (m < inst_.m()) ; ++m ) {
            if (m == -1){
                f << "\"i,1\" -> \"" << inst_.machine(j,0)+1 <<",1\" [ label=\"x(" << j+1 << ",i,1," << inst_.machine(j,0)+1 << ",1)\" ]" << endl;
                continue;
            }
            int m0 = inst_.machine(j,m);
            int mf = (m == inst_.m()-1 ? inst_.m() : inst_.machine(j,m+1));
            int dur = inst_.time(j,m0); // duration time for machine m0 in job j
            for (int t = inst_.est(j,m0); t < inst_.lst(j,m0); t++){
                if (mf < inst_.m()){
                    // arc for another machine
                    f << "\"" << m0+1 << "," << t << "\" -> \"" << mf+1 << "," << t+dur << "\" [ label=\"x(" << j+1 << "," << m0+1 << "," << t << "," << mf+1 << "," << t+dur<<")\" ]" << endl;
                    // arc for same machine (waiting) 
                    //can only be made in the last moment possible
                    if (t == inst_.lst(j,m0)-1) continue;
                    // else
                    f << "\"" << m0+1 << "," << t << "\" -> \"" << m0+1 << "," << t+1 << "\" [ label=\"x(" << j+1 << "," << m0+1 << "," << t << "," << m0+1 << "," << t+1<<")\" ]" << endl;
                } else { // last machine
                    f << "\"" << m0+1 << "," << t << "\" -> \"f," << t+dur << "\" [ label=\"x(" << j+1 << "," << m0+1 << "," << t << ",f," << t+dur<<")\" ]" << endl;
                }
            }
        }
    }


    // string mach0,machf;
    // for ( int j=0 ; (j<inst_.n()) ; ++j ) {
    //     for ( int m0=0 ; (m0<=inst_.m()) ; ++m0 ) {
    //         mach0 = to_string(m0);
    //         if (mach0 == "0") mach0 = "i";
    //         for ( int t0 = 0; t0 < inst_.maxTime()-1; t0++  ) {
    //             for (int mf = 1; mf <= inst_.m()+1; mf++){
    //                 machf = to_string(mf);
    //                 if (machf == to_string(inst_.m()+1)) machf = "f";
    //                 for (int tf=t0+1; tf < inst_.maxTime(); tf++){
    //                     if (m0 == mf && tf > t0+1) continue;
    //                     f << "\"" << m0 << "," << t0+1 << "\" -> " << "\"" << machf << "," << tf << "\" [ label=\"x(" << j+1 << "," << mach0 << "," << t0+1 << "," << machf << "," << tf+1 << ")\" ];"<< endl ; // \" color=\"" << x11_color[j] <<
    //                 }
    //             }
    //         }
    //     }
    // }

    f << "{rank=same;" << " \"i,1\" }" << endl;

    for (int t = 0; t < inst_.maxTime(); t++){
        f << "{rank=same;";
        for (int m = 1; m <= inst_.m()+1; m++){
            string mach;
            if (m == inst_.m()+1) mach = "f";
            else mach = to_string(m);

            if (m > 1 )
                f << ",";
            f << " \"" << mach << "," << t+1 << "\"";
        }
        f << "}" << endl;
    }
    // for (int m = 0; m <= inst_.m()+1; m++){
    //     string mach = to_string(m);
    //     if (mach == "0") mach = "i";
    //     else if (mach == to_string(inst_.m()+1)) mach = "f";
    //     f << "{rank=same;";
    //     for (int t = 0; t < inst_.maxTime(); t++){
    //         f << " \"" << mach << "," << t+1 << "\"";
    //     }
    //     f << "}" << endl;
    // }
    f << "}" << endl;
    f.close();
}