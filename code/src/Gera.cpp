

#include "Gera.hpp"
#include <vector>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>

#include <string.h>

extern "C"
{
    #include "cgraph/lp.h"
    #include "cgraph/strutils.h"
    #include "cgraph/build_cgraph.h"
}

using namespace std;

int Gera::teto(double v){
    if (  fabs(v-floor(v+0.5))<=1e-6 )
        return v;
    else
        return ceil(v);
}


Gera::Gera( const Instance &_inst ) :
    inst_(_inst),
    mip(lp_create()),
    xIdx_(vector< vector< map< int, map< int, map< int, int >>>>> (inst_.n(), vector< map< int, map< int, map< int, int >>>>(inst_.m()+1))),
    enter_Gera(vector<vector<vector<vector<int>>>>(inst_.n(),vector<vector<vector<int>>>(inst_.m()+2+inst_.n(),vector<vector<int>>(inst_.maxTime())))),
    exit_Gera(vector<vector<vector<vector<int>>>>(inst_.n(),vector<vector<vector<int>>>(inst_.m()+1,vector<vector<int>>(inst_.maxTime())))),
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
                exit_Gera[j][0][j].emplace_back(names.size());
                enter_Gera[j][inst_.machine(j,0)+1][0].emplace_back(names.size());
                names.emplace_back( "x("+to_string(j+1)+",i,0,"+to_string(inst_.machine(j,0)+1)+",0)" );
                lb.emplace_back( 0.0 );
                ub.emplace_back( 1 );
                obj.emplace_back( 0 );
                integer.emplace_back( 0 );
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
                    enter_Gera[j][mf+1][t+dur].emplace_back(names.size());
                    exit_Gera[j][m0+1][t].emplace_back(names.size());
                    for (int tp = t; tp < t+dur; tp++){
                        process[j][m0+1][tp].emplace_back(names.size());
                    }
                    names.emplace_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(mf+1)+","+to_string(t+dur)+")" );
                    lb.emplace_back( 0.0 );
                    ub.emplace_back( 1 );
                    obj.emplace_back( 0 );
                    integer.emplace_back( 0 );
                    // arc for same machine (waiting) 
                    // can only be made in the last moment possible
                    if (t == inst_.lst(j,m0)-1) continue;
                    // else
                    xIdx_[j][m0+1][t][m0+1][t+1] = names.size();
                    enter_Gera[j][m0+1][t+1].emplace_back(names.size());
                    exit_Gera[j][m0+1][t].emplace_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.emplace_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(m0+1)+","+to_string(t+1)+")" );
                    lb.emplace_back( 0.0 );
                    ub.emplace_back( 1 );
                    obj.emplace_back( 0 );
                    integer.emplace_back( 0 );
                } else { // conclusion machine f
                    xIdx_[j][m0+1][t][mf+1][t+dur] = names.size();
                    enter_Gera[j][mf+1+j][t+dur].emplace_back(names.size());
                    exit_Gera[j][m0+1][t].emplace_back(names.size());
                    for (int tp = t; tp < t+dur; tp++){
                        process[j][m0+1][tp].emplace_back(names.size());
                    }
                    names.emplace_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+",f,"+to_string(t+dur)+")" );
                    lb.emplace_back( 0.0 );
                    ub.emplace_back( 1 );
                    obj.emplace_back( 0 );
                    integer.emplace_back( 0 );

                    xIdx_[j][m0+1][t][m0+1][t+1] = names.size();
                    enter_Gera[j][m0+1][t+1].emplace_back(names.size());
                    exit_Gera[j][m0+1][t].emplace_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.emplace_back( "x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+","+to_string(m0+1)+","+to_string(t+1)+")" );
                    lb.emplace_back( 0.0 );
                    ub.emplace_back( 1 );
                    obj.emplace_back( 0 );
                    integer.emplace_back( 0 );
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
    
    

    // c var
    cIdx_ = names.size();
    names.emplace_back("C");
    lb.emplace_back( 0.0 );
    ub.emplace_back( inst_.maxTime() );
    obj.emplace_back( 1.0 );
    integer.emplace_back( 0 );
    cout << "Number of variables: " << names.size() << endl;
    clock_t begin = clock();
    lp_add_cols( mip, obj, lb, ub, integer, names );
    clock_t end = clock();
    double time_spent = ((double)end-begin)/((double)CLOCKS_PER_SEC);
    cout << "variáveis criadas. Tempo: " << time_spent << endl;
    // f.open ("exit_Geras.txt");
    // for (int j = 0; j < inst_.n(); j++){
    //     for (int tf=0; tf < inst_.maxTime(); tf++){
    //         for (int mf = 0; mf < inst_.m()+1; mf++){
            
    //             f << "machine " << mf << " time " << tf << endl;
    //             for (int var : exit_Gera[j][mf][tf]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "exit_Geras criado" << endl;

    // f.open ("enter_Geras.txt");
    // for (int j = 0; j < inst_.n(); j++){
    //     for (int t0=0; t0 < inst_.maxTime()-1; t0++){
    //         for (int m0 = 0; m0 <= inst_.m()+1; m0++){
            
    //             f << "machine " << m0 << " time " << t0 << endl;
    //             for (int var : enter_Gera[j][m0][t0]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "enter_Geras criado" << endl;

    // f.open ("processing_machines.txt");
    // for (int t0=0; t0 < inst_.maxTime()-1; t0++){
    //     for (int m0 = 0; m0 <= inst_.m(); m0++){
        
    //         f << "machine " << m0 << " time " << t0 << endl;
    //         for (int j = 0; j < inst_.n(); j++){
    //             for (int var : process[j][m0][t0]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "processing_machines created" << endl;
    // createCompleteGraphDot();
    // cout << "arquivo dot criado" << endl;

    // initial Gera constraint
    vector< string > constr_names;
    vector< int > init_Gera;
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;
        for (int var : exit_Gera[j][0][j]){
            idx.emplace_back( var );
            coef.emplace_back( -1.0 );
        }
        lp_add_row( mip, idx, coef, "init_Gera("+to_string(j)+")", 'E', -1.0 );
        init_Gera.emplace_back(constr_names.size());
        constr_names.emplace_back( "init_Gera("+to_string(j)+")");
    }
    cout << "initial Gera constraints ok" << endl;

    // final Gera constraint
    vector< int > final_Gera;

    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;
        for (int t = 1; t < inst_.maxTime(); t++){
            for (int var : enter_Gera[j][inst_.m()+1+j][t]){
                idx.emplace_back( var );
                coef.emplace_back( 1.0 );
            }
        }

        if (idx.size() != 0){
            lp_add_row( mip, idx, coef, "final_Gera("+to_string(j)+")", 'E', 1.0 );
            final_Gera.emplace_back(constr_names.size());
            constr_names.emplace_back( "final_Gera("+to_string(j)+")");
        }
    }
    cout << "final Gera constraints ok" << endl;
    // Gera constraints
    vector< int > Gera;
    for (int j = 0; j < inst_.n(); j++){
        for ( int t=0 ; (t<inst_.maxTime()) ; ++t ) {
            for ( int m=1 ; (m<=inst_.m()) ; ++m ) {
                vector< int > idx;
                vector< double > coef;
                for (int var : enter_Gera[j][m][t]){
                    idx.emplace_back( var );
                    coef.emplace_back( 1.0 );
                }
                for (int var : exit_Gera[j][m][t]){
                    idx.emplace_back( var );
                    coef.emplace_back( -1.0 );
                }

                if (idx.size() != 0){
                    lp_add_row( mip, idx, coef, "Gera("+to_string(j+1)+","+to_string(m)+","+to_string(t)+")", 'E', 0.0 );
                    Gera.emplace_back(constr_names.size());
                    constr_names.emplace_back("Gera("+to_string(j+1)+","+to_string(m)+","+to_string(t)+")");
                }
            }
        }
    }
    cout << "Gera constraints created" << endl;

    // processing restriction
    vector< int > processing;
    for ( int t=0 ; (t<inst_.maxTime()) ; ++t ) {

        for ( int m=0 ; (m<=inst_.m()) ; ++m ) {
            vector< int > idx;
            vector< double > coef;
            for (int j = 0; j < inst_.n(); j++){
                for( int var : process[j][m][t]){
                    idx.emplace_back( var );
                    coef.emplace_back( 1.0 );

                }
            }
            if (idx.size() != 0){
                lp_add_row( mip, idx, coef, "processing("+to_string(m) + ","+to_string(t) + ")", 'L', 1.0 );        
                processing.emplace_back(constr_names.size());
                constr_names.emplace_back("processing("+to_string(m) + ","+to_string(t) + ")");
            }
        }
        
    }
    cout << "processing constraints created" << endl;

    // restrições fim
    vector< int > fim;
    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        idx.emplace_back( cIdx_ );
        coef.emplace_back( 1.0 );

        for (int t = 1; t < inst_.maxTime(); t++){
            
            for (int var : enter_Gera[j][inst_.m()+1+j][t]){
                idx.emplace_back( var );
                coef.emplace_back( -t );
            }
        }

        lp_add_row( mip, idx, coef, "fim("+to_string(j+1)+")", 'G', 0 );
        fim.emplace_back(constr_names.size());
        constr_names.emplace_back("fim("+to_string(j+1)+")");
    }
    cout << "end constraints created" << endl;

    lp_write_lp( mip, "teste.lp");//inst_.instanceName().c_str() );
    //lp_write_mps( mip, inst_.instanceName().c_str() );

        

}

double Gera::execute(){
    double bnd_continuous = 9;
        int c = 9;
        bool pare = false;
        double limite = 0.00001;
        lp_optimize_as_continuous(mip);
        lp_write_lp(mip,"teste_cb.lp");
        bnd_continuous = lp_obj_value(mip);
        cout << "Continuo = " << bnd_continuous << endl;
        lp_write_lp(mip,"teste_cb.lp");

        lp_as_integer(mip);

        //Callback cb = Callback(mip,inst_,xIdx_,process);
        
        //CGraph *cgraph = build_cgraph_lp(mip);
        //cgraph_print_summary(cgraph, "test_cgraph");
        lp_write_lp(mip,"teste.lp");        
        lp_optimize(mip);
        double bnd_integer = lp_obj_value(mip);

        cout << "Inteiro = " << bnd_integer << endl;
        cout << "Diferença = " << bnd_integer - bnd_continuous << endl;
        //lp_write_lp(mip,"teste_cb.lp");
        lp_write_sol(mip,"teste_cb.sol");
        return bnd_integer - bnd_continuous;
}

Gera::~Gera()
{
    //lp_free( &mip );
}

// void Gera::create_x11(){
//     x11_color.emplace_back("black");
//     x11_color.emplace_back("magenta");
//     x11_color.emplace_back("grey");
//     x11_color.emplace_back("gold");
//     x11_color.emplace_back("blue");
//     x11_color.emplace_back("green");
//     x11_color.emplace_back("yellow");
//     x11_color.emplace_back("red");
//     x11_color.emplace_back("cyan");
//     x11_color.emplace_back("crimson");
//     x11_color.emplace_back("chocolate");
//     x11_color.emplace_back("brown");
//     x11_color.emplace_back("orange");
// }


void Gera::createCompleteGraphDot(){
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

