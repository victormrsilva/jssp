#include "Fernando.hpp"

#include <vector>
#include <string>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <time.h>

#include <string.h>

extern "C"
{
#include "cgraph/lp.h"
#include "cgraph/strutils.h"
#include "cgraph/build_cgraph.h"
}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1, T2> const &pair) const
    {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);

        return h1 ^ h2;
    }
};
using namespace std;

int Fernando::teto(double v)
{
    if (fabs(v - floor(v + 0.5)) <= 1e-6)
        return v;
    else
        return ceil(v);
}

Fernando::Fernando( const Instance &_inst ) : inst_(_inst),
process(vector<vector<vector<vector<int>>>>(inst_.n(), (vector<vector<vector<int>>>(inst_.m(), vector<vector<int>>(inst_.maxTime()+1))))) 
{ // já inicializa a variável privada _inst com o valor passado por referencia
    mip = lp_create();

    // variáveis de decisão
    xIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    eIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    fIdx_ = vector<vector<int>>(inst_.m(),vector<int>(inst_.maxTime()+1));
    enter_flow = vector<vector<vector<vector<int>>>>(inst_.n(), vector<vector<vector<int>>>(inst_.m(), vector<vector<int>>(inst_.maxTime()+2)));

    vector< double > lb; // lower bound
    vector< double > ub; // upper bound
    vector< double > obj; // se é objetivo?
    vector< char > integer; // variável inteira?

    // criação das variáveis x
    for (int i = 0; i < inst_.m(); i++){
        for (int t = 0; t <= inst_.maxTime(); t++){
            fIdx_[i][t] = names.size(); 
            names.push_back("f("+to_string(i+1)+","+to_string(t+1)+")"); // nome dessa variável
            lb.push_back(0.0);
            ub.push_back(1.0);
            obj.push_back(0.0);
            integer.push_back(1);
    
            for (int j = 0; j < inst_.n(); j++){
                int m0 = inst_.machine(j,i);
                if (t >= inst_.est(j,m0) && t <= inst_.lst(j,m0)){
                    xIdx_[m0][j][t] = names.size(); 
                    names.push_back("x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(1);
                    int dur = inst_.time(j, m0);

                    for (int tp = t; tp < t + dur; tp++)
                    {
                        process[j][m0][tp].emplace_back(xIdx_[m0][j][t]);
                    }
                    // int mf = (i == inst_.m() - 1 ? inst_.m() : inst_.machine(j, i + 1));
                    // if (mf != inst_.m()){
                    //     enter_flow[j][mf][t + dur].emplace_back(names.size());
                    // }
                
                    eIdx_[m0][j][t] = names.size(); 
                    names.push_back("e("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(1);
                    // if (mf != inst_.m()){
                    //     enter_flow[j][m0][t + 1].emplace_back(eIdx_[m0][j][t]);
                    // }
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
    //cout << "enter_flows criado" << endl;

    // f.open ("processing_machines.txt");
    // for (int t0=0; t0 < inst_.maxTime()-1; t0++){
    //     for (int m0 = 0; m0 < inst_.m(); m0++){

    //         f << "machine " << m0+1 << " time " << t0 << endl;
    //         for (int j = 0; j < inst_.n(); j++){
    //             for (int var : process[j][m0][t0]){
    //                 f << names[var] << endl;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // cout << "processing_machines created" << endl;

    // adiciona colunas ao solver
    lp_add_cols( mip, obj, lb, ub, integer, names );
    cout << "Number of variables: " << names.size() << endl;

    // restrições
    vector<string> constr_names; // indice das constraints
    //restriction 28
    for (int i = 0; i < inst_.m(); i++){
        vector< int > idx;
        vector< double > coef;

        for (int j = 0; j < inst_.n(); j++){
            if (inst_.machine(j,0) == i) { // se for a primeira máquina
                idx.push_back( xIdx_[i][j][0] );
                coef.push_back( 1.0 );
            }
            // adiciona restrição.
        }
        idx.push_back( fIdx_[i][0] );
        coef.push_back( 1.0 );

        lp_add_row( mip, idx, coef, "inicio_maquina(m"+to_string(i+1)+",t"+to_string(1)+")", 'E', 1 );
        constr_names.emplace_back("inicio_maquina(m"+to_string(i+1)+",t"+to_string(1)+")");
    }
    cout << "restriction inicio_maquina added" << endl;
    // restriction 29
    for (int i = 0; i < inst_.m(); i++){
        for (int t = 1; t <= inst_.maxTime(); t++){
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

            lp_add_row( mip, idx, coef, "fluxo_maquina("+to_string(i+1)+","+to_string(t+1)+")", 'E', 0 );
            constr_names.emplace_back("fluxo_maquina("+to_string(i+1)+","+to_string(t+1)+")");
        }

    }
    cout << "restriction fluxo_maquina added" << endl;

    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        int h = inst_.machine(j,0); // first machine
        cout << h << " " << j << endl;
        idx.push_back( xIdx_[h][j][0] );
        coef.push_back( -1.0 );
        // adiciona restrição.
        idx.push_back( eIdx_[h][j][0] );
        coef.push_back( -1.0 );
        
        lp_add_row( mip, idx, coef, "inicio_espera(m"+to_string(h)+",j"+to_string(j+1)+",t"+to_string(1)+")", 'E', -1 );
        constr_names.emplace_back("inicio_espera(m"+to_string(h)+",j"+to_string(j+1)+",t"+to_string(1)+")");

    }
    cout << "restriction inicio_espera added" << endl;

    for (int i = 0; i < inst_.m(); i++){
        for (int j = 0; j < inst_.n(); j++){
            int h = inst_.machine(j,i); // first machine

            for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
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

                lp_add_row( mip, idx, coef, "fluxo_espera("+to_string(h+1)+","+to_string(j+1)+","+to_string(t+1)+")", 'E', 0 );
                constr_names.emplace_back("fluxo_espera("+to_string(h+1)+","+to_string(j+1)+","+to_string(t+1)+")");
            }
        }
    }
    cout << "restriction fluxo_espera added" << endl;

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
            lp_add_row( mip, idx, coef, "ultimo_tempo("+to_string(h+1)+","+to_string(j+1)+")", 'E', 0.0 );
            constr_names.emplace_back("ultimo_tempo("+to_string(h+1)+","+to_string(j+1)+")");
        }
    }
    cout << "ultimo_tempo constraints created" << endl;

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

        lp_add_row( mip, idx, coef, "fim("+to_string(j+1)+")", 'G', 0 );
        fim.emplace_back(constr_names.size());
        constr_names.emplace_back("fim("+to_string(j+1)+")");
    }
    cout << "makespan constraints created" << endl;

    //lp_write_lp( mip, (inst_.instanceName() + "_machine").c_str() );
    //lp_write_mps( mip, inst_.instanceName().c_str() );
    //lp_optimize_as_continuous(mip);
    if (inst_.execute()){
        optimize();
    }
    // if (inst_.execute()){
    //     lp_optimize( mip );
    //     lp_write_sol(mip, "jssp_Fernando.sol");
    // }
}

void Fernando::cliques(int *idxs,double *coefs)
{
        /* variáveis em conflito
            1 = as variáveis x(j,m(i),t',m(i+1),t'+d) e x(j,m(i),t',m(i),t'+1) com t' > t, ou seja, as variáveis que ainda podem processar e as de espera nos tempos
            2 = as variáveis x(j,m(i-1), t',m(i),t+d) com t'+d > t
            3 = as variáveis de processamento para a máquina i no tempo t
        */

            CGraph *cgraph = cgraph_create(names.size() * 2); // todos os vértices de menos o c
            ofstream file_conflitos("conflitos.txt");

            unordered_set<int> conflitos;
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
                        int idx = xIdx_[m0][j][t];
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
                            // if (tf > inst_.est(j, m0))
                            // {
                            //     file_conflitos << names[xIdx_[j][m0 + 1][tf - 1][m0 + 1][tf]] << " ";
                            //     conflitos.insert(xIdx_[j][m0 + 1][tf - 1][m0 + 1][tf]);
                            // }
                            file_conflitos << names[xIdx_[m0][j][tf]] << " ";
                            conflitos.insert(xIdx_[m0][j][tf]);
                        }
                        file_conflitos << endl;

                        // caso 2
                        file_conflitos << "caso 2: " << endl;
                        if (m != 0)
                        {
                            unordered_set<pair<int,int>,pair_hash> analisar;
                            analisar.emplace(make_pair(m,t));
                            while (!analisar.empty()){
                                pair<int,int> maquina_tempo = *analisar.begin();
                                analisar.erase(analisar.begin());
                                if (maquina_tempo.first == 0) continue; // primeira máquina
                                int m_atual = inst_.machine(j,maquina_tempo.first);
                                int m_anterior = inst_.machine(j, maquina_tempo.first - 1);
                                int t0 = maquina_tempo.second;
                                int dur_anterior = inst_.time(j, m_anterior);
                                if (t0 - dur_anterior >= inst_.time(j, m_anterior))
                                {
                                    for (int tf = t0 - dur_anterior + 1; tf < t0; tf++)
                                    {
                                        cout << "job: " << j+1 << " m_atual: " << m_atual+1 << " t_atual: " << t0 << " m_anterior: " << m_anterior+1 << " dur: " << dur_anterior << " t0: " << tf << endl;
                                        if (tf >= inst_.lst(j, m_anterior))
                                            continue;
                                        cout << names[xIdx_[m_anterior][j][tf]] << endl;
                                        auto i = conflitos.emplace(xIdx_[m_anterior][j][tf]);
                                        if (i.second){ // conseguiu inserir{
                                            file_conflitos << names[xIdx_[m_anterior][j][tf]] << " ";
                                        }
                                        analisar.emplace(make_pair(maquina_tempo.first-1,tf));
                                        
                                    }
                                    cout << analisar.size() << endl;
                                    //getchar();
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
                                conflitos.insert(var);
                            }
                        }
                        file_conflitos << endl;

                        // mostra conflitos
                        file_conflitos << "todos os conflitos: " << endl;
                        file_conflitos << names[idx] << " = ";
                        vector<int> conflicts(conflitos.begin(),conflitos.end());
                        for (int var : conflicts)
                        {
                            file_conflitos << names[var] << " ";
                        }
                        file_conflitos << endl << endl;
                        cgraph_add_node_conflicts(cgraph, idx, &conflicts[0], conflicts.size());
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

void Fernando::optimize(){
    bool continuo = true;
    bool clique = false;
    if (continuo){
        clock_t begin = clock();
        const int nCols = lp_cols(mip);
        int *idxs = new int[nCols];
        double *coefs = new double[nCols];
        lp_optimize_as_continuous(mip);
        int lb = lp_obj_value(mip);
        int ub = inst_.maxTime();
        double bnd = 0;
        double c = ((ub-lb)/2 + lb);
        int iteracoes = 0;
        cout << "lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << bnd << " fabs: " << fabs(lb-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
        //getchar();
        while (fabs(ub-lb) > 1) {
            iteracoes++;
            c = ((ub-lb)/2 + lb);
            bnd = lifting(teto(c),&idxs[0],&coefs[0]);
            //lp_optimize_as_continuous(mip);
            lp_write_lp(mip, "teste_cb.lp");
            lp_write_sol(mip, "solution.sol");
            cout << "antes: lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << bnd << " teto(bnd): " << teto(bnd) << " fabs: " << fabs(c-bnd)  << " fabs: " << fabs(teto(bnd)-c) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
            if (fabs(c-bnd) <= 1e-06) { //lb == floor(bnd)){
                ub = teto(c);
            } else {
                lb = teto(c);
            }
            cout << "depois: lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << bnd << " fabs: " << fabs(c-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
            //getchar();
            // remove colunas de fim
            //            lp_remove_rows(mip,fim);

            //getchar();
        };
        clock_t end = clock();
        cout << "Quantidade de iterações: " << iteracoes << " lb: " << lb << " ub " << ub <<endl;
        double time_spent = ((double)end - begin) / ((double)CLOCKS_PER_SEC);
        cout << "Tempo gasto no lifting: " << time_spent << endl;
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
    string filename = inst_.instanceName()+"_machine";
    lp_write_lp(mip, (filename+".lp").c_str());
    //getchar();
    //lp_optimize(mip);

    //lp_write_lp(mip,"teste_cb.lp");
    lp_write_sol(mip, (filename+".sol").c_str());
}

double Fernando::lifting(int c, int *idxs, double *coefs){
    for (int j = 0; j < inst_.n(); j++){
        vector<int> idx;
        vector<double> coef;
        idx.emplace_back(cIdx_);
        coef.emplace_back(1.0);
        int idxRow = fim[j];
        cout << "idxRow " << idxRow << " name: " << endl;
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
        int h = inst_.machine(j,inst_.m()-1);
        for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
    
            idx.push_back(xIdx_[h][j][t]);
            coef.push_back(-max(t+inst_.time(j,h),c));
        }
        double violado = sol - c;
        cout << endl;
        cout << "C: " << c << " soma: " << sol << " soma - C: " << violado << endl;
        lp_remove_row(mip, fim[j]);
        lp_add_row(mip, idx, coef, "fim(" + to_string(j + 1) + ")", 'G', 0);
    }
    lp_write_lp(mip, "teste_cb.lp");
    
    lp_optimize_as_continuous(mip);
    return lp_obj_value(mip);
}

// void Fernando::optimize(){
//     clock_t begin = clock();
//     bool continuo = true;
//     bool clique = false;
//     if (continuo){
//         const int nCols = lp_cols(mip);
//         int *idxs = new int[nCols];
//         double *coefs = new double[nCols];
//         int lb = 0;
//         int ub = inst_.maxTime();
//         double bnd = 0;
//         int c = 0;
//         int iteracoes = 0;
//         do {
//             iteracoes++;
//             cout << "lb: " << lb << " ub: " << ub << " c: " << c << " iteração: " << iteracoes << endl;
//             lp_optimize_as_continuous(mip);
//             lp_write_lp(mip, "teste_cb.lp");
//             lp_write_sol(mip, "solution.sol");
//             bnd = lp_obj_value(mip);
//             //getchar();
//             // remove colunas de fim
//             //            lp_remove_rows(mip,fim);
//             for (int j = 0; j < inst_.n(); j++){
//                 vector<int> idx;
//                 vector<double> coef;
//                 idx.emplace_back(cIdx_);
//                 coef.emplace_back(1.0);
//                 int idxRow = fim[j];
//                 //cout << "idxRow " << idxRow << endl;
//                 const int nElements = lp_row(mip, idxRow, idxs, coefs);
//                 double sol = 0;

//                 for (int i = 0; i < nElements; i++){
//                     char *nome = lp_varName(mip, idxs[i]);
//                     if (strcmp(nome, "C") == 0)
//                         continue;
//                     double x = lp_xIdx(mip, idxs[i]);
//                     //lp_row_name(mip,idxs[i], nome);
//                     int c2 = max(-1 * (int)coefs[i], c);
//                     //cout << coefs[i] << " " << c2 << " " <<  idxs[i] << " " << nome << " " << x << " " << c2*x <<endl;
//                     sol += c2 * x;
//                 }

//                 int h = inst_.machine(j,inst_.m()-1);
//                 for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
            
//                     idx.push_back(xIdx_[h][j][t]);
//                     coef.push_back(-max(t+inst_.time(j,h),c));
//                 }   

//                 double violado = sol - c;
//                 cout << endl;
//                 cout << "C: " << c << " soma: " << sol << " soma - C: " << violado << endl;
//                 lp_remove_row(mip, fim[j]);
//                 lp_add_row(mip, idx, coef, "fim(" + to_string(j + 1) + ")", 'G', 0);
//             }
            
//             if (c == teto(bnd)){
//                 ub = teto(bnd);
//             } else {
//                 lb = teto(bnd);
//             }
//             c = teto((ub-lb)/2 + lb);
//             cout << "lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << teto(bnd) << endl;
//             lp_write_lp(mip, "teste_cb.lp");
//             getchar();
//         } while (ub != lb);
//         cout << "Quantidade de iterações: " << iteracoes << endl;
//         clock_t end = clock();
//         double time_spent = ((double)end - begin) / ((double)CLOCKS_PER_SEC);
//         cout << "Tempo gasto no lifting: " << time_spent << endl;
//         getchar();
//         if (clique)
//         {
//             cliques(idxs,coefs);
//         }
//     }

//     lp_as_integer(mip);

//     //Callback cb = Callback(mip,inst_,xIdx_,process);

//     //CGraph *cgraph = build_cgraph_lp(mip);
//     //cgraph_print_summary(cgraph, "test_cgraph");
//     lp_write_lp(mip, "teste.lp");
//     //getchar();
//     lp_optimize(mip);

//     //lp_write_lp(mip,"teste_cb.lp");
//     lp_write_sol(mip, "teste_cb.sol");
// }

Fernando::~Fernando()
{
    lp_free( &mip );
}

