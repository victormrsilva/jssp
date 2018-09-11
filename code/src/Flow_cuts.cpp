

#include "Flow_cuts.hpp"
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

int Flow::teto(double v)
{
    if (fabs(v - floor(v + 0.5)) <= 1e-6)
        return v;
    else
        return ceil(v);
}

Flow::Flow(const Instance &_inst) : inst_(_inst),
                                    mip(lp_create()),
                                    xIdx_(vector<vector<map<int, map<int, map<int, int>>>>>(inst_.n(), vector<map<int, map<int, map<int, int>>>>(inst_.m() + 1))),
                                    enter_flow(vector<vector<vector<vector<int>>>>(inst_.n(), vector<vector<vector<int>>>(inst_.m() + 2 + inst_.n(), vector<vector<int>>(inst_.maxTime()+1)))),
                                    exit_flow(vector<vector<vector<vector<int>>>>(inst_.n(), vector<vector<vector<int>>>(inst_.m() + 1, vector<vector<int>>(inst_.maxTime()+1)))),
                                    process(vector<vector<vector<vector<int>>>>(inst_.n(), (vector<vector<vector<int>>>(inst_.m() + 1, vector<vector<int>>(inst_.maxTime()+1)))))
{
    vector<double> lb;
    vector<double> ub;
    vector<double> obj;
    vector<char> integer;
    
    // creating x vars
    for (int j = 0; (j < inst_.n()); ++j)
    {
        for (int m = -1; (m < inst_.m()); ++m)
        {
            if (m == -1)
            { // máquina inicial
                xIdx_[j][m + 1][0][inst_.machine(j, 0) + 1][0] = names.size();
                exit_flow[j][0][j].emplace_back(names.size());
                enter_flow[j][inst_.machine(j, 0) + 1][0].emplace_back(names.size());
                names.emplace_back("x(" + to_string(j + 1) + ",i,0," + to_string(inst_.machine(j, 0) + 1) + ",0)");
                lb.emplace_back(0.0);
                ub.emplace_back(1);
                obj.emplace_back(0);
                integer.emplace_back(1);
                continue;
            }
            int m0 = inst_.machine(j, m);
            int mf = (m == inst_.m() - 1 ? inst_.m() : inst_.machine(j, m + 1));
            int dur = inst_.time(j, m0); // duration time for machine m0 in job j
            for (int t = inst_.est(j, m0); t <= inst_.lst(j, m0); t++)
            {
                if (mf < inst_.m())
                {
                    // arc for another machine
                    xIdx_[j][m0 + 1][t][mf + 1][t + dur] = names.size();
                    //cout << j << " " << m0+1 << " " << t << " " << mf+1 << " " << t+dur << endl;
                    enter_flow[j][mf + 1][t + dur].emplace_back(names.size());
                    exit_flow[j][m0 + 1][t].emplace_back(names.size());
                    for (int tp = t; tp < t + dur; tp++)
                    {
                        process[j][m0 + 1][tp].emplace_back(names.size());
                    }
                    //names.emplace_back("x(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t) + "," + to_string(mf + 1) + "," + to_string(t + dur) + ")");
                    names.emplace_back("x(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t) + ")");
                    lb.emplace_back(0.0);
                    ub.emplace_back(1);
                    obj.emplace_back(0);
                    integer.emplace_back(1);
                    // arc for same machine (waiting)
                    // can only be made in the last moment possible
                    if (t == inst_.lst(j, m0))
                        continue;
                    // else
                    xIdx_[j][m0 + 1][t][m0 + 1][t + 1] = names.size();
                    enter_flow[j][m0 + 1][t + 1].emplace_back(names.size());
                    exit_flow[j][m0 + 1][t].emplace_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.emplace_back("e(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t) + ")");
                    // names.emplace_back("x(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t) + "," + to_string(m0 + 1) + "," + to_string(t + 1) + ")");
                    lb.emplace_back(0.0);
                    ub.emplace_back(1);
                    obj.emplace_back(0);
                    integer.emplace_back(1);
                }
                else
                { // conclusion machine f
                    xIdx_[j][m0 + 1][t][mf + 1][t + dur] = names.size();
                    enter_flow[j][mf + 1 + j][t + dur].emplace_back(names.size());
                    exit_flow[j][m0 + 1][t].emplace_back(names.size());
                    for (int tp = t; tp < t + dur; tp++)
                    {
                        process[j][m0 + 1][tp].emplace_back(names.size());
                    }
                    names.emplace_back("x(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t)  + ")");
                    lb.emplace_back(0.0);
                    ub.emplace_back(1);
                    obj.emplace_back(0);
                    integer.emplace_back(1);

                    xIdx_[j][m0 + 1][t][m0 + 1][t + 1] = names.size();
                    enter_flow[j][m0 + 1][t + 1].emplace_back(names.size());
                    exit_flow[j][m0 + 1][t].emplace_back(names.size());
                    //cout << j << " " << m0+1 << " " << t << " " << m0+1 << " " << t+1 << endl;
                    names.emplace_back("e(" + to_string(j + 1) + "," + to_string(m0 + 1) + "," + to_string(t) + ")");
                    lb.emplace_back(0.0);
                    ub.emplace_back(1);
                    obj.emplace_back(0);
                    integer.emplace_back(1);
                }
            }
        }
    }

    ofstream f;
    f.open("variables.txt");
    for (string name : names)
    {
        f << name << endl;
    }
    f.close();

    // c var
    cIdx_ = names.size();
    names.emplace_back("C");
    lb.emplace_back(0.0);
    ub.emplace_back(inst_.maxTime()+1);
    obj.emplace_back(1.0);
    integer.emplace_back(1);
    cout << "Number of variables: " << names.size() << endl;
    clock_t begin = clock();
    lp_add_cols(mip, obj, lb, ub, integer, names);
    clock_t end = clock();
    double time_spent = ((double)end - begin) / ((double)CLOCKS_PER_SEC);
    cout << "variáveis criadas. Tempo: " << time_spent << endl;
    // f.open ("exit_flows.txt");
    // for (int j = 0; j < inst_.n(); j++){
    //     for (int tf=0; tf < inst_.maxTime(); tf++){
    //         for (int mf = 0; mf < inst_.m()+1; mf++){

    //             f << "machine " << mf << " time " << tf << endl;
    //             for (int var : exit_flow[j][mf][tf]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "exit_flows criado" << endl;

    // f.open ("enter_flows.txt");
    // for (int j = 0; j < inst_.n(); j++){
    //     for (int t0=0; t0 < inst_.maxTime()-1; t0++){
    //         for (int m0 = 0; m0 <= inst_.m()+1; m0++){

    //             f << "machine " << m0 << " time " << t0 << endl;
    //             for (int var : enter_flow[j][m0][t0]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "enter_flows criado" << endl;

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

    // initial flow constraint
    vector<string> constr_names;
    vector<int> init_flow;
    for (int j = 0; j < inst_.n(); j++)
    {
        vector<int> idx;
        vector<double> coef;
        for (int var : exit_flow[j][0][j])
        {
            idx.emplace_back(var);
            coef.emplace_back(1.0);
        }
        lp_add_row(mip, idx, coef, "init_flow(" + to_string(j) + ")", 'E', 1.0);
        init_flow.emplace_back(constr_names.size());
        constr_names.emplace_back("init_flow(" + to_string(j) + ")");
    }
    cout << "initial flow constraints ok" << endl;

    // final flow constraint
    vector<int> final_flow;

    for (int j = 0; j < inst_.n(); j++)
    {
        vector<int> idx;
        vector<double> coef;
        for (int t = 1; t <= inst_.maxTime(); t++)
        {
            for (int var : enter_flow[j][inst_.m() + 1 + j][t])
            {
                idx.emplace_back(var);
                coef.emplace_back(1.0);
            }
        }

        if (idx.size() != 0)
        {
            lp_add_row(mip, idx, coef, "final_flow(" + to_string(j) + ")", 'E', 1.0);
            final_flow.emplace_back(constr_names.size());
            constr_names.emplace_back("final_flow(" + to_string(j) + ")");
        }
    }
    cout << "final flow constraints ok" << endl;
    // flow constraints
    vector<int> flow;
    for (int j = 0; j < inst_.n(); j++)
    {
        for (int t = 0; (t <= inst_.maxTime()); ++t)
        {
            for (int m = 1; (m <= inst_.m()); ++m)
            {
                vector<int> idx;
                vector<double> coef;
                for (int var : enter_flow[j][m][t])
                {
                    idx.emplace_back(var);
                    coef.emplace_back(1.0);
                }
                for (int var : exit_flow[j][m][t])
                {
                    idx.emplace_back(var);
                    coef.emplace_back(-1.0);
                }

                if (idx.size() != 0)
                {
                    lp_add_row(mip, idx, coef, "flow(" + to_string(j + 1) + "," + to_string(m) + "," + to_string(t) + ")", 'E', 0.0);
                    flow.emplace_back(constr_names.size());
                    constr_names.emplace_back("flow(" + to_string(j + 1) + "," + to_string(m) + "," + to_string(t) + ")");
                }
            }
        }
    }
    cout << "flow constraints created" << endl;

    // processing restriction
    vector<int> processing;
    for (int t = 0; (t <= inst_.maxTime()); ++t)
    {

        for (int m = 0; (m <= inst_.m()); ++m)
        {
            vector<int> idx;
            vector<double> coef;
            for (int j = 0; j < inst_.n(); j++)
            {
                for (int var : process[j][m][t])
                {
                    idx.emplace_back(var);
                    coef.emplace_back(1.0);
                }
            }
            if (idx.size() != 0)
            {
                lp_add_row(mip, idx, coef, "processing(" + to_string(m) + "," + to_string(t) + ")", 'L', 1.0);
                processing.emplace_back(constr_names.size());
                constr_names.emplace_back("processing(" + to_string(m) + "," + to_string(t) + ")");
            }
        }
    }
    cout << "processing constraints created" << endl;

    // restrições fim

    for (int j = 0; j < inst_.n(); j++){
        vector<int> idx;
        vector<double> coef;

        idx.emplace_back(cIdx_);
        coef.emplace_back(1.0);

        for (int t = 1; t <= inst_.maxTime(); t++){

            for (int var : enter_flow[j][inst_.m() + 1 + j][t]){
                idx.emplace_back(var);
                coef.emplace_back(-t);
            }
        }

        lp_add_row(mip, idx, coef, "fim(" + to_string(j + 1) + ")", 'G', 0);
        fim.emplace_back(constr_names.size());
        constr_names.emplace_back("fim(" + to_string(j + 1) + ")");
    }
    cout << "end constraints created" << endl;

    lp_write_lp(mip, (inst_.instanceName() + "_flow").c_str()); //inst_.instanceName().c_str() );
    //lp_optimize(mip);
    //lp_write_sol(mip, "jssp_Flow.sol");

    //lp_write_mps( mip, inst_.instanceName().c_str() );
//    lp_optimize_as_continuous(mip);
    if (inst_.execute()){
        optimize();
    }
}

void Flow::cliques(int *idxs,double *coefs)
{
        /* variáveis em conflito
            1 = as variáveis x(j,m(i),t',m(i+1),t'+d) e x(j,m(i),t',m(i),t'+1) com t' > t, ou seja, as variáveis que ainda podem processar e as de espera nos tempos
            2 = as variáveis x(j,m(i-1), t',m(i),t+d) com t'+d > t
            3 = as variáveis de processamento para a máquina i no tempo t
        */

            CGraph *cgraph = cgraph_create(names.size() * 2); // todos os vértices de menos o c
            ofstream file_conflitos("conflitos.txt");

            vector<int> conflitos;
            //cgraph_add_node_conflicts(cgraph,cIdx_,&conflitos[0],conflitos.size());
            vector<int> indices_conflitos;
            for (int m = 0; m < inst_.m(); m++)
            {
                for (int j = 0; j < inst_.n(); j++)
                {
                    int m0 = inst_.machine(j, m);
                    int mf = (m == inst_.m() - 1 ? inst_.m() : inst_.machine(j, m + 1));
                    int dur = inst_.time(j, m0); // duration time for machine m0 in job j
                    for (int t = inst_.est(j, m0); t <= inst_.lst(j, m0); t++)
                    {
                        //cout << j << " " << m0 << " " << t << " " << mf << " " << t+dur << endl;
                        int idx = xIdx_[j][m0 + 1][t][mf + 1][t + dur];
                        // cout << idx << endl;
                        indices_conflitos.emplace_back(idx);
                        conflitos.clear();
                        file_conflitos << "variavel: " << idx << " " << names[idx] << endl;
                        // caso 1
                        file_conflitos << "caso 1: " << endl;
                        for (int tf = inst_.est(j, m0); tf < inst_.lst(j, m0); tf++)
                        {
                            if (t == tf)
                                continue;
                            if (tf > inst_.est(j, m0))
                            {
                                file_conflitos << names[xIdx_[j][m0 + 1][tf - 1][m0 + 1][tf]] << " ";
                                conflitos.emplace_back(xIdx_[j][m0 + 1][tf - 1][m0 + 1][tf]);
                            }
                            file_conflitos << names[xIdx_[j][m0 + 1][tf][mf + 1][tf + dur]] << " ";
                            conflitos.emplace_back(xIdx_[j][m0 + 1][tf][mf + 1][tf + dur]);
                        }
                        file_conflitos << endl;

                        // caso 2
                        file_conflitos << "caso 2: " << endl;
                        if (m != 0)
                        {
                            vector<pair<int,int>> analisar;
                            analisar.emplace_back(make_pair(m,t));
                            while (!analisar.empty()){
                                pair<int,int> maquina_tempo = analisar.back();
                                analisar.pop_back();
                                if (maquina_tempo.first == 0) continue; // primeira máquina
                                int m_anterior = inst_.machine(j, maquina_tempo.first - 1);
                                int t0 = maquina_tempo.second;
                                int dur_anterior = inst_.time(j, m_anterior);
                                if (t0 - dur_anterior >= inst_.time(j, m_anterior))
                                {
                                    for (int tf = t0 - dur_anterior + 1; tf < t0; tf++)
                                    {
                                        cout << "job: " << j+1 << " m_atual: " << maquina_tempo.first+1 << " t_atual: " << t0 << " m_anterior: " << m_anterior+1 << " dur: " << dur_anterior << " t0: " << tf << endl;
                                        cout << names[xIdx_[j][m_anterior + 1][tf][m0 + 1][tf + dur_anterior]] << endl;
                                        if (tf >= inst_.lst(j, m_anterior))
                                            continue;

                                        file_conflitos << names[xIdx_[j][m_anterior + 1][tf][m0 + 1][tf + dur_anterior]] << " ";
                                        conflitos.emplace_back(xIdx_[j][m_anterior + 1][tf][m0 + 1][tf + dur_anterior]);
                                        analisar.emplace_back(make_pair(maquina_tempo.first-1,tf));
                                        
                                    }
                                    cout << analisar.size() << endl;
                                    getchar();
                                }
                            }
                        }
                        // caso 3
                        file_conflitos << endl;
                        file_conflitos << "caso 3: " << endl;
                        for (int j_ = 0; j_ < inst_.n(); j_++)
                        {
                            for (int var : process[j_][m0 + 1][t])
                            {
                                if (var == idx)
                                    continue;
                                file_conflitos << names[var] << " ";
                                conflitos.emplace_back(var);
                            }
                        }
                        file_conflitos << endl
                                        << endl;

                        // mostra conflitos
                        file_conflitos << "todos os conflitos: " << endl;
                        file_conflitos << names[idx] << " = ";
                        for (int var : conflitos)
                        {
                            file_conflitos << names[var] << " ";
                        }
                        file_conflitos << endl;

                        cgraph_add_node_conflicts(cgraph, idx, &conflitos[0], conflitos.size());
                        // idx = xIdx_[j][m0+1][t][m0+1][t+1];
                        // conflitos.clear();
                        // cgraph_add_node_conflicts(cgraph,idx,&conflitos[0],conflitos.size());
                        //cout << "conflitos adicionados no cgraph" << endl;
                        //getchar();
                    }
                }
            }
            file_conflitos.close();
            cgraph_save(cgraph, "cgraph.txt");
            cout << indices_conflitos.size() << endl;

            double *x = lp_x(mip);
            double *rc = lp_reduced_cost(mip);

            for (unsigned int i = 0; i < names.size(); i++)
            {
                if (x[i] != 0 || rc[i] != 0)
                    cout << names[i] << " " << i << " " << x[i] << " " << rc[i] << endl;
            }

            cout << endl;
            vector<double> x_conflitos = vector<double>(names.size() * 2);
            vector<double> rc_conflitos = vector<double>(names.size() * 2);

            for (unsigned int i = 0; i < names.size(); i++)
            {
                x_conflitos[i] = x[i];
                x_conflitos[i + names.size()] = 1 - x[i];
                rc_conflitos[i] = rc[i];
                rc_conflitos[i + names.size()] = (-1) * rc[i];
            }

            CliqueSeparation *clique_sep = clq_sep_create(cgraph);
            cout << "clique_sep ok" << endl;
            clq_sep_set_verbose(clique_sep, 'T');
            clq_sep_set_rc(clique_sep, &rc_conflitos[0]); //&rc_conflitos[0]);
            cout << "clique_sep_set_rc ok" << endl;
            getchar();
            clq_sep_separate(clique_sep, &x_conflitos[0]); //&x_conflitos[0]);
            cout << "clique_separate ok" << endl;
            getchar();
            const CliqueSet *cliques = clq_sep_get_cliques(clique_sep);
            clq_set_print(cliques);
            getchar();
            int qtd_cliques = clq_set_number_of_cliques(cliques);
            ofstream file_cliques("cliques.txt");
            file_cliques << "qtd de cliques: " << qtd_cliques << endl;

            for (int i = 0; i < qtd_cliques; i++)
            {
                const IntSet *clq = clq_set_get_clique(cliques, i);
                file_cliques << "clique " << i << " tamanho " << clq->size << endl;
                file_cliques << "elementos: ";
                for (int j = 0; j < clq->size; j++)
                {
                    file_cliques << names[clq->elements[j]] << "[" << clq->elements[j] << "], " ;
                }
                file_cliques << endl;
            }
            file_cliques.close();
            delete []idxs;
            delete []coefs;
}

Flow::~Flow()
{
    //lp_free( &mip );
}

// void Flow::create_x11(){
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

void Flow::createCompleteGraphDot()
{
    ofstream f;
    f.open("complete.dot");

    f << "digraph complete {" << endl;
    //f << "ratio = \"auto\" ;" << endl;
    f << "rankdir=LR;" << endl
      << endl;
    for (int m = 0; m <= inst_.m() + 1; m++)
    {
        if (m == 0)
        {
            f << "\"i,1\" [shape=box];" << endl;
            continue;
        }
        string mach;
        if (mach == to_string(inst_.m() + 1))
            mach = "f";
        else
            mach = to_string(m);
        for (int t = 0; t <= inst_.maxTime(); t++)
        {
            f << "\"" << mach << "," << t + 1 << "\" [shape=box];" << endl;
        }
    }

    for (int j = 0; (j < inst_.n()); ++j)
    {
        for (int m = -1; (m < inst_.m()); ++m)
        {
            if (m == -1)
            {
                f << "\"i,1\" -> \"" << inst_.machine(j, 0) + 1 << ",1\" [ label=\"x(" << j + 1 << ",i,1," << inst_.machine(j, 0) + 1 << ",1)\" ]" << endl;
                continue;
            }
            int m0 = inst_.machine(j, m);
            int mf = (m == inst_.m() - 1 ? inst_.m() : inst_.machine(j, m + 1));
            int dur = inst_.time(j, m0); // duration time for machine m0 in job j
            for (int t = inst_.est(j, m0); t <= inst_.lst(j, m0); t++)
            {
                if (mf < inst_.m())
                {
                    // arc for another machine
                    f << "\"" << m0 + 1 << "," << t << "\" -> \"" << mf + 1 << "," << t + dur << "\" [ label=\"x(" << j + 1 << "," << m0 + 1 << "," << t << "," << mf + 1 << "," << t + dur << ")\" ]" << endl;
                    // arc for same machine (waiting)
                    //can only be made in the last moment possible
                    if (t == inst_.lst(j, m0))
                        continue;
                    // else
                    f << "\"" << m0 + 1 << "," << t << "\" -> \"" << m0 + 1 << "," << t + 1 << "\" [ label=\"x(" << j + 1 << "," << m0 + 1 << "," << t << "," << m0 + 1 << "," << t + 1 << ")\" ]" << endl;
                }
                else
                { // last machine
                    f << "\"" << m0 + 1 << "," << t << "\" -> \"f," << t + dur << "\" [ label=\"x(" << j + 1 << "," << m0 + 1 << "," << t << ",f," << t + dur << ")\" ]" << endl;
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

    f << "{rank=same;"
      << " \"i,1\" }" << endl;

    for (int t = 0; t <= inst_.maxTime(); t++)
    {
        f << "{rank=same;";
        for (int m = 1; m <= inst_.m() + 1; m++)
        {
            string mach;
            if (m == inst_.m() + 1)
                mach = "f";
            else
                mach = to_string(m);

            if (m > 1)
                f << ",";
            f << " \"" << mach << "," << t + 1 << "\"";
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

void Flow::optimize(){
    bool continuo = true;
    bool clique = true;
    if (continuo){
        const int nCols = lp_cols(mip);
        int *idxs = new int[nCols];
        double *coefs = new double[nCols];
        int bnd_anterior = 0;
        double bnd = 9;
        int c = 9;
        do {
            bnd_anterior = c;
            lp_optimize_as_continuous(mip);
            lp_write_lp(mip, "teste_cb.lp");
            lp_write_sol(mip, "solution.sol");
            bnd = lp_obj_value(mip);
            c = teto(bnd);
            // remove colunas de fim
            //            lp_remove_rows(mip,fim);
            for (int j = 0; j < inst_.n(); j++){
                vector<int> idx;
                vector<double> coef;
                idx.emplace_back(cIdx_);
                coef.emplace_back(1.0);
                int idxRow = fim[j];
                //cout << "idxRow " << idxRow << endl;
                const int nElements = lp_row(mip, idxRow, idxs, coefs);
                double sol = 0;

                for (int i = 0; i < nElements; i++){
                    char *nome = lp_varName(mip, idxs[i]);
                    if (strcmp(nome, "C") == 0)
                        continue;
                    double x = lp_xIdx(mip, idxs[i]);
                    //lp_row_name(mip,idxs[i], nome);
                    int c2 = max(-1 * (int)coefs[i], c);
                    //cout << coefs[i] << " " << c2 << " " <<  idxs[i] << " " << nome << " " << x << " " << c2*x <<endl;
                    sol += c2 * x;
                }

                for (int t = 1; t <= inst_.maxTime(); t++){
                    for (int var : enter_flow[j][inst_.m() + 1 + j][t]){
                        //int coeficiente = max(t,c);
                        idx.emplace_back(var);
                        coef.emplace_back(-max(t, c));
                    }
                }
                double violado = sol - c;
                cout << endl
                        << "C: " << c << " soma: " << sol << " soma - C: " << violado << endl;
                lp_remove_row(mip, fim[j]);
                lp_add_row(mip, idx, coef, "fim(" + to_string(j + 1) + ")", 'G', 0);
            }
            lp_write_lp(mip, "teste_cb.lp");
            //getchar();
        } while (bnd_anterior != c);
        //getchar();

        if (clique)
        {
            cliques(idxs,coefs);
        }
    }

    lp_as_integer(mip);

    //Callback cb = Callback(mip,inst_,xIdx_,process);

    //CGraph *cgraph = build_cgraph_lp(mip);
    //cgraph_print_summary(cgraph, "test_cgraph");
    lp_write_lp(mip, "teste.lp");
    //getchar();
    lp_optimize(mip);

    //lp_write_lp(mip,"teste_cb.lp");
    lp_write_sol(mip, "teste_cb.sol");
}