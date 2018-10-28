#include "Fernando.hpp"

#include <vector>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <unordered_set>

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

double Fernando::teto(double v)
{
    //cout << v << endl;
    if ( fabs(v-(ceil(v)))<=1e-10 )
        return v;
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
            names.push_back("f("+to_string(i+1)+","+to_string(t)+")"); // nome dessa variável
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
    ub.push_back( inst_.maxTime() );
    obj.push_back( 1.0 );
    integer.push_back( 1 );

    ofstream f;
    f.open ("variables_pack.txt");
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

    variables_pack = vector<vector<int>>(names.size());
    cout << variables_pack.size() << endl;
    
    for (int m0 = 0; m0 < inst_.m(); m0++){
        for (int j = 0; j < inst_.n(); j++){
            for (int t=inst_.est(j,m0); t <= inst_.lst(j,m0); t++){
                int dur = inst_.time(j,m0);
                int var = xIdx_[m0][j][t];
                cout << names[var] << " ";
                for (int t0 = t - inst_.time(j,m0) + 1; t0 < t+inst_.minimumTime(j); t0++){
                    if (t0 < inst_.est(j,m0) || t0 > inst_.lst(j,m0) || t0 == t) continue;
                    cout << "[" << j+1 << "," << m0+1 << "," << t0 << "," << inst_.est(j,m0) << "," << inst_.lst(j,m0) << "] ";
                    cout << names[xIdx_[m0][j][t0]] << " ";
                    variables_pack[var].emplace_back(xIdx_[m0][j][t0]);
                }
                
                for (int j_aux = 0; j_aux < inst_.n(); j_aux++){
                    if (j_aux == j) continue;
                    int dur_aux = inst_.time(j_aux,m0);
                    int qtd_movimentos = inst_.minimumTime(j) - dur_aux ;
                    if (qtd_movimentos >= 0){
                        while (qtd_movimentos >= 0){
                            int tempo = t+qtd_movimentos;
                            if (tempo >= inst_.est(j_aux,m0) && tempo <= inst_.lst(j_aux,m0)){
                                cout << "[" << j_aux+1 << "," << m0+1 << "," << tempo << "," << inst_.est(j_aux,m0) << "," << inst_.lst(j_aux,m0) << "] ";
                                cout << names[xIdx_[m0][j_aux][tempo]] << " ";
                                variables_pack[var].emplace_back(xIdx_[m0][j_aux][tempo]);
                            }
                            qtd_movimentos--;
                        }
                    }
                    else {
                        while (qtd_movimentos <= 0){
                            int tempo = t+qtd_movimentos;
                            if (tempo >= inst_.est(j_aux,m0) && tempo <= inst_.lst(j_aux,m0)){
                                cout << "[" << j_aux+1 << "," << m0+1 << "," << tempo << "," << inst_.est(j_aux,m0) << "," << inst_.lst(j_aux,m0) << "] ";
                                cout << names[xIdx_[m0][j_aux][tempo]] << " ";
                                variables_pack[var].emplace_back(xIdx_[m0][j_aux][tempo]);
                            }
                            qtd_movimentos++;
                        }
                    }
                }
                cout << endl;
            }
        }
    }

    f.open("cuts.txt");
    for (int i = 0; i < variables_pack.size(); i++){
        f << names[i] << " = ";
        for (int v : variables_pack[i]){
            f << names[v] << " ";
        }
        f << endl;
    }
    f.close();
    getchar();
    f.open ("processing_machines.txt");
    for (int t0=0; t0 < inst_.maxTime()-1; t0++){
        for (int m0 = 0; m0 < inst_.m(); m0++){

            f << "machine " << m0+1 << " time " << t0 << endl;
            for (int j = 0; j < inst_.n(); j++){
                for (int var : process[j][m0][t0]){
                    f << names[var] << endl;
                }
            }
        }
    }
    f.close();
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

            lp_add_row( mip, idx, coef, "fluxo_maquina("+to_string(i+1)+","+to_string(t)+")", 'E', 0 );
            constr_names.emplace_back("fluxo_maquina("+to_string(i+1)+","+to_string(t)+")");
        }

    }
    cout << "restriction fluxo_maquina added" << endl;

    for (int j = 0; j < inst_.n(); j++){
        vector< int > idx;
        vector< double > coef;

        int h = inst_.machine(j,0); // first machine
        cout << h << " " << j << endl;
        idx.push_back( xIdx_[h][j][0] );
        coef.push_back( 1.0 );
        // adiciona restrição.
        idx.push_back( eIdx_[h][j][0] );
        coef.push_back( 1.0 );
        
        lp_add_row( mip, idx, coef, "inicio_espera(m"+to_string(h)+",j"+to_string(j+1)+",t"+to_string(1)+")", 'E', 1 );
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

                lp_add_row( mip, idx, coef, "fluxo_espera("+to_string(h+1)+","+to_string(j+1)+","+to_string(t)+")", 'E', 0 );
                constr_names.emplace_back("fluxo_espera("+to_string(h+1)+","+to_string(j+1)+","+to_string(t)+")");
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

    // for (int i = 0; i < inst_.m(); i++){
    //     for (int j = 0; j < inst_.n(); j++){
    //         for (int t = inst_.est(j,i); t <= inst_.lst(j,i); t++){
    //             //cout << j << " " << i << " " << t << endl;
    //             vector< int > idx;
    //             vector< double > coef;
    //             idx.push_back(fIdx_[i][t]);
    //             coef.push_back(1.0);
    //             for (int tf = t; tf <= inst_.lst(j,i); tf++){
    //                 idx.push_back(xIdx_[i][j][tf]);
    //                 coef.push_back(1.0);
    //             }
    //             lp_add_row( mip, idx, coef, "left_align("+to_string(j+1)+","+to_string(i+1)+","+to_string(t)+")", 'L', 1.0 );
    //             constr_names.emplace_back("left_align("+to_string(j+1)+","+to_string(i+1)+","+to_string(t)+")");
    //         }
    //     }
    // }

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

    lp_write_lp( mip, (inst_.instanceName() + "_machine_original.lp").c_str() );
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

int Fernando::manual_cuts(){
    double *x = lp_x(mip);
    ofstream f;
    f.open("cuts_manual.txt");
    int qtd = 0;
    for (int i = 0; i < variables_pack.size(); i++){
        f << names[i] << "(" << x[i] << ") ";
        double soma = x[i]; // inicia a soma para ver se há violação
        for (int var : variables_pack[i]){
            f << names[var] << "(" << x[var] << ") ";
            soma = soma + x[var];
        }
        f << soma << endl;
        
        if (soma > 1.001){ // houve violação então vamos adicionar os cortes
            cout << soma << endl;
            vector< int > idx;
            vector< double > coef;
            idx.emplace_back(i);
            coef.emplace_back(1.0);
            for (int var : variables_pack[i]){
                idx.emplace_back(var);
                coef.emplace_back(1.0);
            }
            if (idx.size() > 1){
                qtd_manual_cuts++;
                qtd++;
                lp_add_row( mip, idx, coef, "manual_cut("+to_string(qtd_manual_cuts)+")", 'L', 1 );
            }
        }
    }
    f.close();
    string filename = inst_.instanceName()+"_cortes";
    lp_write_lp(mip, (filename+".lp").c_str());
    getchar();
    return qtd;
}

int Fernando::cliques(int *idxs,double *coefs)
{
    clock_t begin = clock();
    double *x = lp_x(mip);
    double *rc = lp_reduced_cost(mip);

    // for (unsigned int i = 0; i < names.size(); i++)
    // {
    //     if (x[i] != 0 || rc[i] != 0)
    //         cout << names[i] << " " << i << " " << x[i] << " " << rc[i] << endl;
    // }

    // cout << endl;
    vector<double> x_conflitos = vector<double>(names.size() * 2);
    vector<double> rc_conflitos = vector<double>(names.size() * 2);
    double delta = 1e-6;
    for (unsigned int i = 0; i < names.size(); i++)
    {
        x_conflitos[i] = (x[i] == 0 ? delta : x[i]);
        x_conflitos[i + names.size()] = 1 - (x[i] == 0 ? delta : x[i]);
        rc_conflitos[i] = rc[i];
        rc_conflitos[i + names.size()] = (-1) * rc[i];
    }

    CliqueSeparation *clique_sep = clq_sep_create(cgraph);
    cout << "clique_sep ok" << endl;
    clq_sep_set_verbose(clique_sep, 'T');
    clq_sep_set_rc(clique_sep, &rc_conflitos[0]); //&rc_conflitos[0]);
    cout << "clique_sep_set_rc ok" << endl;
    //getchar();
    clq_sep_separate(clique_sep, &x_conflitos[0]); //&x_conflitos[0]);
    cout << "clique_separate ok" << endl;
    //getchar();
    const CliqueSet *cliques = clq_sep_get_cliques(clique_sep);
    //clq_set_print(cliques);
    //getchar();
    int qtd_cliques = clq_set_number_of_cliques(cliques);
    ofstream file_cliques("cliques.txt");
    file_cliques << "qtd de cliques: " << qtd_cliques << endl;

    for (int i = 0; i < qtd_cliques; i++)
    {
        vector< int > idx;
        vector< double > coef;

        const IntSet *clq = clq_set_get_clique(cliques, i);
        file_cliques << "clique " << i << " tamanho " << clq->size << endl;
        file_cliques << "elementos: ";
        for (int j = 0; j < clq->size; j++)
        {
            idx.push_back( clq->elements[j] );
            coef.push_back( 1.0 );
            file_cliques << names[clq->elements[j]] << " + " ;
        }
        file_cliques << " <= 1 " << endl;
        lp_add_row( mip, idx, coef, "cortes("+to_string(qtd_cortes)+")", 'L', 1 );
        qtd_cortes++;
    }
    file_cliques.close();
    string filename = inst_.instanceName()+"_cortes";
    lp_write_lp(mip, (filename+".lp").c_str());
    clock_t end = clock();
    cout << "cuts added: " << qtd_cliques << " time for adding on lp: " << (double) (end-begin)/CLOCKS_PER_SEC << endl;
    cout << "file cliques.txt created. LP " << filename << " created for this iteration. Press enter to continue" << endl;
    return qtd_cliques;
}
void Fernando::cgraph_creation()
{
    /* variáveis em conflito
    1 = as variáveis x(j,m(i),t',m(i+1),t'+d) e x(j,m(i),t',m(i),t'+1) com t' > t, ou seja, as variáveis que ainda podem processar e as de espera nos tempos
    2 = as variáveis x(j,m(i-1), t',m(i),t+d) com t'+d > t
    3 = as variáveis de processamento para a máquina i no tempo t*/

    clock_t begin = clock();
    cgraph = cgraph_create(names.size() * 2); // todos os vértices de menos o c
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

                        // caso 2
                        file_conflitos << endl << "caso 2: " << endl;
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
                                        //cout << "job: " << j+1 << " m_atual: " << m_atual+1 << " t_atual: " << t0 << " m_anterior: " << m_anterior+1 << " dur: " << dur_anterior << " t0: " << tf << endl;
                                        if (tf >= inst_.lst(j, m_anterior))
                                            continue;
                                        //cout << names[xIdx_[m_anterior][j][tf]] << endl;
                                        auto i = conflitos.emplace(xIdx_[m_anterior][j][tf]);
                                        if (i.second){ // conseguiu inserir{
                                            file_conflitos << names[xIdx_[m_anterior][j][tf]] << " ";
                                        }
                                        analisar.emplace(make_pair(maquina_tempo.first-1,tf));
                                        
                                    }
                                    //cout << analisar.size() << endl;
                                    //getchar();
                                }
                            }
                        }
                        
                        // caso 3
                        file_conflitos << endl << "caso 3: " << endl;
                        for (int j_ = 0; j_ < inst_.n(); j_++)
                        {
                            for (int var : process[j_][m0][t])
                            {
                                if (var == idx)
                                    continue;
                                //file_conflitos << names[var] << " ";
                                conflitos.insert(var);
                            }
                        }

                        // caso 4
                        file_conflitos << endl << "caso 4: " << endl;
                        if (m != 0)
                        {
                            unordered_set<pair<int,int>,pair_hash> analisar;
                            //cout << "analisar j: " << j+1 << " maquina: " << inst_.machine(j,m-1)+1 << " tempo: " << t-inst_.time(j, inst_.machine(j,m-1)) << " maquina atual: " << inst_.machine(j,m)+1 << " tempo atual: " << t << endl;
                            analisar.emplace(make_pair(m-1,t-inst_.time(j, inst_.machine(j,m-1))));
                            while (!analisar.empty()){
                                pair<int,int> maquina_tempo = *analisar.begin();
                                analisar.erase(analisar.begin());
                                if (maquina_tempo.first == -1) continue; // primeira máquina
                                int t0 = maquina_tempo.second;
                                int m0 = inst_.machine(j,maquina_tempo.first);
                                int mf = inst_.machine(j,maquina_tempo.first-1);
                                //cout << "j:" << j+1 << " m: " << m << " t: " << t << " m0: " << m0+1 << " mf: " << mf+1 << " t0: " << t0 << " tf: " << t0-inst_.time(j,mf) << endl;
                                //cout << maquina_tempo.first-1 << " " << t0-inst_.time(j,mf) << endl;
                                for (int tf = t0 + 1 ; tf <= inst_.lst(j,m0); tf++)
                                {
                                    //cout << names[xIdx_[m0][j][tf]] << endl;
                                    auto i = conflitos.emplace(xIdx_[m0][j][tf]);
                                    if (i.second){ // conseguiu inserir{
                                        file_conflitos << names[xIdx_[m0][j][tf]] << " ";
                                    }
                                }
                                
                                analisar.emplace(make_pair(maquina_tempo.first-1,t0-inst_.time(j,mf)));
                                //cout << analisar.size() << endl;
                                //getchar();
                            }
                        }
                        // mostra conflitos
                        file_conflitos << endl << "todos os conflitos: " << endl;
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
    clock_t end = clock();
    cout << "cgraph creation time: " << (double) (end-begin)/CLOCKS_PER_SEC << endl;
    //cgraph_save(cgraph, "cgraph.txt");
    //cout << indices_conflitos.size() << endl;
}

void Fernando::optimize(){

    if (continuo){
        const int nCols = lp_cols(mip);
        int *idxs = new int[nCols];
        double *coefs = new double[nCols];
        if (binario){
            lifting_binario(idxs,coefs);
        } else {
            lifting_linear(idxs,coefs);
        }
        //getchar();
        delete []idxs;
        delete []coefs;
    }
    getchar();
    lp_as_integer(mip);

    //Callback cb = Callback(mip,inst_,xIdx_,process);

    //CGraph *cgraph = build_cgraph_lp(mip);
    //cgraph_print_summary(cgraph, "test_cgraph");
    string filename = inst_.instanceName()+"_machine";
    lp_write_lp(mip, (filename+".lp").c_str());
    //getchar();
    lp_optimize(mip);

    //lp_write_lp(mip,"teste_cb.lp");
    lp_write_sol(mip, (filename+".sol").c_str());
}

double Fernando::lifting(double c, int *idxs, double *coefs){
    string filename = inst_.instanceName()+"_lifting";
    for (int j = 0; j < inst_.n(); j++){
        vector<int> idx;
        vector<double> coef;
        idx.emplace_back(cIdx_);
        coef.emplace_back(1.0);
        int idxRow = lp_get_constr_by_name(mip,("fim("+to_string(j+1)+")").c_str());
        //cout << "idxRow " << idxRow << " name: " << endl;
        const int nElements = lp_row(mip, idxRow, idxs, coefs);
        double sol = 0;
        double aux = 0;
        for (int i = 0; i < nElements; i++){
            char *nome = lp_varName(mip, idxs[i]);
            if (strcmp(nome, "C") == 0)
                continue;
            double x = lp_xIdx(mip, idxs[i]);
            //lp_row_name(mip,idxs[i], nome);
            if ((-1 * coefs[i]) > c )
                aux = -1 * coefs[i];
            else 
                aux = c;
            //cout << coefs[i] << " " << c2 << " " <<  idxs[i] << " " << nome << " " << x << " " << c2*x <<endl;
            sol += aux * x;
        }
        int h = inst_.machine(j,inst_.m()-1);
        for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
    
            idx.push_back(xIdx_[h][j][t]);
            if ((fabs(t+inst_.time(j,h)) - c ) > 1e-6)
                aux = t+inst_.time(j,h);
            else 
                aux = c;
            coef.push_back(-aux);
        }
        double violado = sol - c;
        cout << endl;
        cout << "C: " << c << " soma: " << sol << " soma - C: " << violado << endl;
        lp_remove_row(mip, idxRow);
        lp_add_row(mip, idx, coef, "fim(" + to_string(j + 1) + ")", 'G', 0);
        lp_write_lp(mip, "teste_cb.lp");
    }
    lp_write_lp(mip, (filename+".lp").c_str());
    
    lp_optimize_as_continuous(mip);
    lp_write_sol(mip,(filename+".sol").c_str());
    cout << "Lifting solved for c = " << c << ". Files " << filename+".lp" << " e " << filename+".sol" << " successfully created." << endl;
    //getchar();
    return lp_obj_value(mip);
}

void Fernando::lifting_binario(int *idxs, double *coefs){
    if (clique){
        cgraph_creation();
    }
    ofstream saida(inst_.instanceName()+"_solution_binario_machine.csv");
    saida << "iteracao;lb;ub;c;cortes;valor;tempo" << endl;
    clock_t begin = clock();
    lp_optimize_as_continuous(mip);
    double lb = teto(lp_obj_value(mip));
    double ub = inst_.maxTime();
    saida << "0;0;" << ub << ";0;" << lp_obj_value(mip) <<";"<<lp_solution_time(mip)<<endl;
    double bnd = 0;
    double c = ((ub-lb)/2 + lb);
    int iteracoes = 0;
    cout << "lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << bnd << " fabs: " << fabs(lb-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
    //getchar();
    while (fabs(ub-lb) > 1e-6) {
        iteracoes++;
        c = ((ub-lb)/2 + lb);
        saida << iteracoes << ";" << lb << ";" << ub << ";" <<teto(c) << ";";
        if (clique)
        {
            manual_cuts();
            //saida << cliques(idxs,coefs) << ";";
        }
        
        bnd = lifting(teto(c),&idxs[0],&coefs[0]);
        saida << bnd << ";" << lp_solution_time(mip) << endl;
        //lp_optimize_as_continuous(mip);
        //lp_write_lp(mip, "teste_cb.lp");
        //lp_write_sol(mip, "solution.sol");
        cout << "antes: lb: " << lb << " ub: " << ub << " c: " << c << " bnd: " << bnd << " teto(bnd): " << teto(bnd) << " teto(c): " << teto(c)   << endl; //<< " floor(bnd):" << floor(bnd) << endl;
        if (fabs(teto(c) - bnd) < 1e-6) { //
            ub = teto(c);
        } else {
            lb = teto(c);
        }
        cout << "depois: lb: " << lb << " ub: " << ub << " teto(c): " << teto(c) << " bnd: " << bnd << " fabs: " << fabs(c-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
        
        if (fabs(ub - lb) <= 1){
            cout << "antes: lb: " << lb << " ub: " << ub << " teto(c): " << teto(c) << " bnd: " << bnd << " fabs: " << fabs(c-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
            c = ((ub-lb)/2 + lb);
            c = teto(c);
            bnd = lifting(teto(c),&idxs[0],&coefs[0]);
            if (fabs(c - bnd) < 1e-6) { //
                lb = c;
            } else {
                ub = c;
            }
            cout << "depois: lb: " << lb << " ub: " << ub << " teto(c): " << teto(c) << " bnd: " << bnd << " fabs: " << fabs(c-bnd) << endl; //<< " floor(bnd):" << floor(bnd) << endl;
        }

        //getchar();
        // remove colunas de fim
        //            lp_remove_rows(mip,fim);

        //getchar();
    };
    
    clock_t end = clock();
    cout << "Quantidade de iterações binario: " << iteracoes << " lb: " << lb << " ub " << ub <<endl;
    saida << iteracoes << ";" << lb << ";" << ub << ";" << teto(c) << ";-;-" << endl;
    saida.close();
    double time_spent = ((double)end - begin) / ((double)CLOCKS_PER_SEC);
    cout << "Tempo gasto no lifting binario: " << time_spent << endl;
    string filename = inst_.instanceName()+"_machine_lift_bin";
    lp_write_lp(mip, (filename+".lp").c_str());
    lp_write_sol(mip, (filename+".sol").c_str());

}

void Fernando::lifting_linear(int *idxs, double *coefs){
    if (clique){
        cgraph_creation();
    }
    clock_t begin = clock();
    int iteracoes = 0;
    lp_optimize_as_continuous(mip);
    double bnd = lp_obj_value(mip);
    double bnd_anterior = 0;
    double c = 0;
    ofstream saida(inst_.instanceName()+"_solution_linear_machine.csv");
    saida << "iteracao;bnd;c;cortes;tempo" << endl;
    saida << "0;" << bnd << ";0;0;"<<lp_solution_time(mip)<<endl;
    while (fabs(bnd - teto(bnd_anterior)) > 1e-06) {
        iteracoes++;
        int cortes = 0;
        if (clique)
        {
            cortes = manual_cuts();
            cortes = cortes + cliques(idxs,coefs);
        }
        
        bnd_anterior = bnd;
        c = teto(bnd);
        bnd = lifting(c,idxs,coefs);
        
        //lp_write_lp(mip, "teste_cb.lp");
        //lp_write_sol(mip, "solution.sol");
        saida << iteracoes << ";" << bnd << ";" << c << ";" << cortes << ";" << lp_solution_time(mip)<<endl;
        cout <<  " c: " << c << " bnd_anterior: " << bnd_anterior << " bnd: " << bnd << endl;
        getchar();
    } 
    clock_t end = clock();
    cout << "Quantidade de iterações linear: " << iteracoes  <<endl;
    double time_spent = ((double)end - begin) / ((double)CLOCKS_PER_SEC);
    cout << "Tempo gasto no lifting linear: " << time_spent << endl;
    string filename = inst_.instanceName()+"_machine_lift_lin";
    lp_write_lp(mip, (filename+".lp").c_str());
    lp_write_sol(mip, (filename+".sol").c_str());
}


Fernando::~Fernando()
{
    lp_free( &mip );
}

