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



using std::cout;
using std::vector;
using std::ofstream;
using std::endl;
using std::unordered_set;
using std::string;
using std::to_string;
using std::pair;

bool compare(int i, int j) { return i < j;}

double Fernando::teto(double v)
{
    //cout << v << endl;
    if ( fabs(v-(ceil(v)))<=1e-10 )
        return v;
    return ceil(v);
}

Fernando::Fernando( Instance &_inst ) : inst_(_inst),
process(vector<vector<vector<vector<int>>>>(inst_.n(), (vector<vector<vector<int>>>(inst_.m(), vector<vector<int>>(inst_.maxTime()+1))))) 
{ // já inicializa a variável privada _inst com o valor passado por referencia
    // reduz_lst_kondili(5);

    buildProblem();

    buildCliqueCuts();

    lp_write_lp( mip, (inst_.instanceName() + "_machine_original.lp").c_str() );

    bool criterio_parada = false;
    double bnd_anterior = 0;
    lp_optimize(mip);
    double bnd = lp_obj_value(mip);
    int iteracao = 1;
    ofstream f;
    f.open(inst_.instanceName() + "_continuos_fenchel_3.csv");
    f << "iteracao;bnd;cliques;fenchel;valor_lifting" << endl;
    f << "0;" << bnd <<";0;0;0" << endl;
    while (!criterio_parada){
        cout << "iteracao " << iteracao << endl;
        int cliqueIter = 0;
        int clique = 0;
        int fenchel = 0;
        int fenchelIter = 0;
        lp_set_print_messages(mip,0);
        do {            
            clique = manual_cuts();
            fenchel = executeFenchel();
            cliqueIter += clique;
            fenchelIter += fenchel;
            cout << "clique: " << clique << " fenchel: " << fenchel << endl;
            lp_write_lp(mip,"lp_after_cuts.lp");
            lp_optimize(mip);
        } while ((clique+fenchel) > 0);
        lp_set_print_messages(mip,1);
        cout << "qtd cortes: clique: " << cliqueIter << " fenchel: " << fenchelIter << endl;
        bnd_anterior = bnd;
        // cuts
        // int clique = manual_cuts();
        // int fenchel = executeFenchel();
        // getchar();
        int c = teto(bnd);
        // getchar();
        bnd = lifting(c); // calculate lifting. the optimization is made within this function
        f << iteracao << ";" << bnd << ";" << cliqueIter << ";" << fenchelIter << ";" << c << endl;
        cout << bnd << " " << bnd_anterior << " " << teto(bnd_anterior)<< " " << cliqueIter << " " << fenchelIter << endl;
        // getchar();
        if (fabs(bnd - teto(bnd_anterior)) < 1e-06){
            cout << bnd << " " << teto(bnd_anterior) << endl;
            criterio_parada = true;
            f.close();
        }
        iteracao++;
    }
    lp_write_lp(mip,(inst_.instanceName() + "_continuos_cuts.lp").c_str());

    //lp_write_mps( mip, inst_.instanceName().c_str() );
    //lp_optimize_as_continuous(mip);
    // if (inst_.execute()){
    //     optimize();
    // }
    // // if (inst_.execute()){
    //     lp_optimize( mip );
    //     lp_write_sol(mip, "jssp_Fernando.sol");
    // }
}

int Fernando::manual_cuts(){
    double *x = lp_x(mip);
    // ofstream f;
    // f.open("cuts_manual.txt");
    int qtd = 0;
    for (vector<int> aux : variables_pack){
        //f << names[i] << "(" << x[i] << ") ";
        double soma = 0; // inicia a soma para ver se há violação
        for (int var : aux){
            // f << names[var] << "(" << x[var] << ") ";
            soma = soma + x[var];
        }
        // f << soma << endl;
        
        if (soma > 1.001){ // houve violação então vamos adicionar os cortes
            //cout << soma << endl;
            vector< int > idx;
            vector< double > coef;
            for (int var : aux){
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
    // f.close();
    // string filename = inst_.instanceName()+"_cortes";
    // lp_write_lp(mip, (filename+".lp").c_str());
    //getchar();
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
    cout << "clique_sep";
    CliqueSeparation *clique_sep = clq_sep_create(cgraph);
    cout << " ok" << endl;
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
    // ofstream file_cliques("cliques.txt");
    // file_cliques << "qtd de cliques: " << qtd_cliques << endl;

    for (int i = 0; i < qtd_cliques; i++)
    {
        vector< int > idx;
        vector< double > coef;

        const IntSet *clq = clq_set_get_clique(cliques, i);
        // file_cliques << "clique " << i << " tamanho " << clq->size << endl;
        // file_cliques << "elementos: ";
        for (int j = 0; j < clq->size; j++)
        {
            idx.push_back( clq->elements[j] );
            coef.push_back( 1.0 );
            // file_cliques << names[clq->elements[j]] << " + " ;
        }
        // file_cliques << " <= 1 " << endl;
        lp_add_row( mip, idx, coef, "cortes("+to_string(totalCliques)+")", 'L', 1 );
        totalCliques++;
    }
    // delete []x;
    // delete []rc;
    // file_cliques.close();
    string filename = inst_.instanceName()+"_cortes";
    lp_write_lp(mip, (filename+".lp").c_str());
    clock_t end = clock();
    cout << "cuts added: " << qtd_cliques << " time for adding on lp: " << (double) (end-begin)/CLOCKS_PER_SEC << endl;
    // cout << "file cliques.txt created. LP " << filename << " created for this iteration. Press enter to continue" << endl;
    return qtd_cliques;
}

int Fernando::oddHoles(){
    int numCols = lp_cols(mip);
    CutPool *cutPool = cut_pool_create(numCols);

    int oddCuts = lp_generate_odd_hole_cuts(mip, cgraph, cutPool);
    ofstream f("oddholes.txt");
    f << "Oddholes encontrados: " << oddCuts << endl;
    for(int j = 0; j < oddCuts; j++)
    {
        string name = "oddHole("+to_string(j)+")";
        const Cut *cut = cut_pool_get_cut(cutPool, j);
        for (int i = 0; i < cut_size(cut); i++){
            f << names[cut_get_idxs(cut)[i]] << " ";
        }
        f << " < " << cut_get_rhs(cut) << endl;
        lp_add_row(mip, cut_size(cut), cut_get_idxs(cut), cut_get_coefs(cut), name.c_str(), 'L', cut_get_rhs(cut));
        f.close();
    }

    cut_pool_free(&cutPool);
    return oddCuts;
}
void Fernando::cgraph_creation()
{
    /* variáveis em conflito
    1 = as variáveis x(j,m(i),t',m(i+1),t'+d) e x(j,m(i),t',m(i),t'+1) com t' > t, ou seja, as variáveis que ainda podem processar e as de espera nos tempos
    2 = as variáveis x(j,m(i-1), t') com t' entre t-t(j,m(i-1)) e est(j,m(i-1))
    3 = as variáveis de processamento para a máquina i no tempo t
    4 = as variáveis x(j,m(i+1),t'), com t' < t+d t*/

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
                            //cout << "analisar j: " << j+1 << " maquina: " << inst_.machine(j,m-1)+1 << " tempo: " << t-inst_.time(j, inst_.machine(j,m-1)) << " maquina atual: " << inst_.machine(j,m)+1 << " tempo atual: " << t << endl;
                            analisar.emplace(std::make_pair(m-1,t-inst_.time(j, inst_.machine(j,m-1))));
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
                                
                                analisar.emplace(std::make_pair(maquina_tempo.first-1,t0-inst_.time(j,mf)));
                                //cout << analisar.size() << endl;
                                //getchar();
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
                                file_conflitos << names[var] << " ";
                                conflitos.insert(var);
                            }
                        }

                        // caso 4
                        file_conflitos << endl << "caso 4: " << endl;


                        for (int k = m+1; k < inst_.m(); k++){
                            int m0 = inst_.machine(j,m);
                            int mf = inst_.machine(j,k);
                            int m_anterior = inst_.machine(j,k-1);
                            //cout << "j" << j << " m0: " << m0+1 << " t0: " << t0 << " mf: " << mf+1 << " k: " << k << " distance: " << inst_.distance(j,m0,mf) << " est(mf): " << inst_.est(j,mf) << endl;
                            for (int tf = inst_.est(j,mf); tf < t+inst_.distance(j,m0,mf); tf++){
                                //cout << m0+1 << " " << t0 << " " << mf+1 << " " << tf << " " << t0+inst_.time(j,m_anterior) << endl;
                                int var = xIdx_[mf][j][tf];
                                conflitos.insert(var);
                                file_conflitos << names[var] << " ";
                                //cout << names[var] << endl;
                            }
                        }
                                    //getchar();
                        
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
    //getchar();
    //cgraph_save(cgraph, "cgraph.txt");
    //cout << indices_conflitos.size() << endl;
}

void Fernando::optimize(){

    // if (continuo){
    //     const int nCols = lp_cols(mip);
    //     int *idxs = new int[nCols];
    //     double *coefs = new double[nCols];
    //     if (binario){
    //         lifting_binario(idxs,coefs);
    //     } else {
    //         lifting_linear(idxs,coefs);
    //     }
    //     //getchar();
    //     delete []idxs;
    //     delete []coefs;
    // }
    // //getchar();
    // string filename = inst_.instanceName()+"_machine_cuts_lift_continuos";
    // lp_optimize_as_continuous(mip);
    // lp_write_lp(mip, (filename+".lp").c_str());
    // lp_write_sol(mip, (filename+".sol").c_str());

    //Callback cb = Callback(mip,inst_,xIdx_,process);

    //CGraph *cgraph = build_cgraph_lp(mip);
    //cgraph_print_summary(cgraph, "test_cgraph");
    // filename = inst_.instanceName()+"_machine_integer";
    // lp_write_lp(mip, (filename+".lp").c_str());
    //getchar();
    //lp_optimize(mip);

    //lp_write_lp(mip,"teste_cb.lp");
    //lp_write_sol(mip, (filename+".sol").c_str());
}

// calculate the new lifting and optmize as continuous
double Fernando::lifting(double c){
    int *idxs = new int[lp_cols(mip)];
    double *coefs = new double[lp_cols(mip)];

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
    
    lp_optimize(mip);
    lp_write_sol(mip,(filename+".sol").c_str());
    cout << "Lifting solved for c = " << c << ". Files " << filename+".lp" << " e " << filename+".sol" << " successfully created." << endl;
    // getchar();
    delete []idxs;
    delete []coefs;
    return lp_obj_value(mip);
}

void Fernando::lifting_binario(){
    // if (clique){
    //     cgraph_creation();
    // }
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
        int cortes = 0;
        c = ((ub-lb)/2 + lb);
        
        if (clique)
        {
            cortes = manual_cuts();
            //saida << cliques(idxs,coefs) << ";";
        }
        
        bnd = lifting(teto(c));
        saida << iteracoes << ";" << lb << ";" << ub << ";" <<teto(c) << ";" << cortes <<"," << bnd << ";" << lp_solution_time(mip) << endl;
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
            bnd = lifting(teto(c));
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

void Fernando::lifting_linear(){
    // if (clique){
    //     cgraph_creation();
    // }
    clock_t begin = clock();
    int iteracoes = 0;
    lp_optimize_as_continuous(mip);
    double bnd = lp_obj_value(mip);
    double bnd_anterior = 0;
    double c = 0;
    ofstream saida(inst_.instanceName()+"_solution_linear_machine_cortes.csv");
    saida << "iteracao;bnd;c;cortes;tempo" << endl;

    int cortes = 0;
    saida << 0 << ";" << bnd << ";" << c << ";" << cortes << ";" << lp_solution_time(mip)<<endl;

    while (fabs(bnd - teto(bnd_anterior)) > 1e-06) {
        iteracoes++;
        cortes = 0;
        if (clique)
        {
            cortes = manual_cuts();
            //cortes = oddHoles();
            //cout << "oddHoles: " << cortes << endl;
            //cortes = cortes + cliques(idxs,coefs);
        }
        
        bnd_anterior = bnd;
        c = teto(bnd);
        bnd = lifting(c);
        
        //lp_write_lp(mip, "teste_cb.lp");
        //lp_write_sol(mip, "solution.sol");
        saida << iteracoes << ";" << bnd << ";" << c << ";" << cortes << ";" << lp_solution_time(mip)<<endl;
        cout <<  " c: " << c << " bnd_anterior: " << bnd_anterior << " bnd: " << bnd << endl;
        //getchar();
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

void Fernando::combinacao(int job, unsigned int tam, vector<int> &vec, vector<vector<int> > &combinacoes){
    if (vec.size() == tam) return;
    for (int j = 0; j < inst_.n(); j++){
        if (j != job && find(vec.begin(), vec.end(), j) == vec.end()){
            vec.push_back(j);
            if (vec.size() < tam){
                combinacao(job, tam,vec,combinacoes);
            }
                
            if (vec.size() == tam){
                if (combinacoes.size()==0){
                    // cout << "tam 0: " << job << " ";
                    // for (int i : vec){
                    //     cout << i << " ";
                    // }
                    // cout << endl;
                    combinacoes.emplace_back(vec);
                }
                else {
                    bool existe = false;
                    for(vector<int> v : combinacoes){
                        if (existe){
                            break;
                        }
                        if (isSubset(v,vec)){
                            existe = true;
                            break;
                        }
                    }
                    if (!existe){
                        // cout << "tam " << combinacoes.size() << ": " << job << " ";
                        // for (int i : vec){
                        //     cout << i << " ";
                        // }
                        // cout << endl;
                        combinacoes.emplace_back(vec);
                        // getchar();
                    }
                }
            }
            vec.pop_back();
        }
    }
}

void Fernando::reduz_lst_kondili(int k_max){
    int k = 1;
    if (k_max > inst_.n()){
        k_max = inst_.n();
    }

    while (k < k_max){
        cout << k << " " << k_max << endl;
        bool mudou = false;
        vector<vector<int>> combinacoes;
        for (int i = 0; i < inst_.n(); i++){
            vector<vector<int>> aux;
            vector<int> vec;
            //vec.push_back(i);
            combinacao(i, k ,vec,aux);
            for (vector<int> v : aux){
                v.insert(v.begin(),i);
                combinacoes.emplace_back(v);
            }
        }

        cout << "combinacoes: " << combinacoes.size() << endl;

        for(vector<int> aux : combinacoes){
            string nome = "";
            for (unsigned int i = 0; i < aux.size(); i++){
                nome += to_string(aux[i]) + "-";
            }

            cout << nome << " ";
            LinearProgram *mip_teste = lp_create();
            lp_set_print_messages(mip_teste,0);
            vector< vector< vector< int > > > xIdx_teste;

            xIdx_teste = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime())));

            vector< string > names_teste; // nome das variáveis
            vector< double > lb; // lower bound
            vector< double > ub; // upper bound
            vector< double > obj; // se é objetivo?
            vector< char > integer; // variável inteira?
            

            // criação das variáveis x
            for (int i = 0; i < inst_.m(); i++){
                for (unsigned int j1= 0; j1 < aux.size(); j1++){
                    int j = aux[j1];
                    for (int t = inst_.est(j,i); t <= inst_.lst(j,i); t++){
                        xIdx_teste[i][j][t] = names_teste.size(); // número do índice da variável (vai de 1 até n*m)
                        names_teste.push_back("x("+to_string(j+1)+","+to_string(i+1)+","+to_string(t)+")"); // nome dessa variável
                        //cout << names_teste[xIdx_teste[i][j][t]] << " ";
                        lb.push_back(0.0);
                        ub.push_back(DBL_MAX);
                        obj.push_back(0.0);
                        integer.push_back(1);
                    }
                }
            }
            //cout << endl;
            // c var
            int cIdx_teste = names_teste.size();
            names_teste.push_back("C");
            lb.push_back( 0.0 );
            ub.push_back( DBL_MAX );
            obj.push_back( -1.0 );
            integer.push_back( 1 );

            // adiciona colunas ao solver
            lp_add_cols( mip_teste, obj, lb, ub, integer, names_teste );
            // cout << "var ok" << endl;
            // cout << "Number of variables: " << names_teste.size() << endl;


            // restrição de tempo
            for (int i = 0; i < inst_.m(); i++){
                for (unsigned int j1= 0; j1 < aux.size(); j1++){
                    int j = aux[j1];
                    vector< int > idx;
                    vector< double > coef;
                    for (int t = inst_.est(j,i); t <= inst_.lst(j,i); t++){
                        idx.push_back( xIdx_teste[i][j][t] );
                        coef.push_back( 1.0 );
                        // adiciona restrição.
                    }
                    lp_add_row( mip_teste, idx, coef, "time("+to_string(i+1)+","+to_string(j+1)+")", 'E', 1 );
                }
            }
            // cout << "time restriction added" << endl;
            for (int i = 0; i < inst_.m(); i++){
                for (int t = 0; t < inst_.maxTime(); t++){
                    vector< int > idx;
                    vector< double > coef;

                    for (unsigned int j1= 0; j1 < aux.size(); j1++){
                        int j = aux[j1];
                        //bool mudou = false;
                        for (int t_aux = 0; t_aux < inst_.time(j,i); t_aux++){
                            if ((t - t_aux) >= inst_.est(j,i) && (t-t_aux) <= inst_.lst(j,i)){
                                idx.push_back( xIdx_teste[i][j][t - t_aux] );
                                coef.push_back( 1 );
                            }
                        }
                    }

                    lp_add_row( mip_teste, idx, coef, "processing("+to_string(i+1)+","+to_string(t)+")", 'L', 1 );
                }
            }
            // cout << "processing restriction added" << endl;
            for (unsigned int j1= 0; j1 < aux.size(); j1++){
                int j = aux[j1];
                for (int i = 1; i < inst_.m(); i++){
                    int h = inst_.machine(j,i);
                    vector< int > idx;
                    vector< double > coef;

                    for (int t = inst_.est(j,h); t <= inst_.lst(j,h); t++){
                        int h1 = inst_.machine(j,i-1);
                        int t_aux = t - inst_.time(j,h1);
                        if (t_aux >= inst_.est(j,h1) && t_aux <= inst_.lst(j,h1)){
                            idx.push_back( xIdx_teste[h1][j][t_aux] );
                            coef.push_back( t  );
                        }
                        idx.push_back( xIdx_teste[h][j][t] );
                        coef.push_back( -1 * t );
                    }
                    lp_add_row( mip_teste, idx, coef, "ord("+to_string(h+1)+","+to_string(j+1)+")", 'L', 0 );
                    
                }
            }
            // cout << "ord restriction added" << endl;
            // makespan
            int j = aux[0];
            vector< int > idx;
            vector< double > coef;
            idx.push_back( cIdx_teste );
            coef.push_back( 1 );

            for (int t = inst_.est(j,inst_.machine(j,0)); t <= inst_.lst(j,inst_.machine(j,0)); t++){
                //cout << names_teste[xIdx_teste[inst_.machine(j,0)][j][t]] << " ";
                idx.push_back( xIdx_teste[inst_.machine(j,0)][j][t] );
                coef.push_back( -1 * (t+inst_.time(j,inst_.machine(j,0))) );
                // adiciona restrição.
            }
            //cout << endl;
            lp_add_row( mip_teste, idx, coef, "makespan("+to_string(inst_.machine(j,0)+1)+","+to_string(j+1)+")", 'E', 0 );
            
            // cout << "makespan constraints created" << endl;
            lp_set_print_messages(mip_teste,0);
            // lp_write_lp( mip_teste, (inst_.instanceName() + "_"+nome+".lp").c_str() );
            lp_optimize( mip_teste );
            // lp_write_sol(mip_teste, (inst_.instanceName() + "_"+ nome+".sol").c_str() );
            int maquina = inst_.machine(aux[0],0);
            int lst = -1*lp_obj_value(mip_teste) - inst_.time(aux[0],maquina);
            cout << "lst antigo: " << inst_.lst(aux[0],maquina) << " lst novo: " << lst << endl;
            //getchar();
            if (lst < inst_.lst(aux[0],maquina)){
                inst_.setLst(aux[0],maquina,lst);
                cout << "lst modificado: " << inst_.lst(aux[0],maquina) << endl;
                mudou = true;
                inst_.saveCmpl("jssp.cdat");
                getchar();
            }

            lp_free( &mip_teste );
            //getchar();
        }
        if (mudou){
            k = 2;
        } else {
            k++;
        }
    }
}

void Fernando::buildProblem(){
    mip = lp_create();

    // variáveis de decisão
    xIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    eIdx_ = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1)));
    fIdx_ = vector<vector<int>>(inst_.m(),vector<int>(inst_.maxTime()+1));

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
            integer.push_back(0);
    
            for (int j = 0; j < inst_.n(); j++){
                int m0 = inst_.machine(j,i);
                if (t >= inst_.est(j,m0) && t <= inst_.lst(j,m0)){
                    xIdx_[m0][j][t] = names.size(); 
                    names.push_back("x("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(0);
                    int dur = inst_.time(j, m0);

                    for (int tp = t; tp < t + dur; tp++) {
                        process[j][m0][tp].emplace_back(xIdx_[m0][j][t]);
                    }
                
                    eIdx_[m0][j][t] = names.size(); 
                    names.push_back("e("+to_string(j+1)+","+to_string(m0+1)+","+to_string(t)+")"); // nome dessa variável
                    lb.push_back(0.0);
                    ub.push_back(1.0);
                    obj.push_back(0.0);
                    integer.push_back(0);
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
    integer.push_back( 0 );

    ofstream f;

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
        //cout << h << " " << j << endl;
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

    // last possible time
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
}

void Fernando::buildCliqueCuts(){
    cout << variables_pack.size() << endl;
    
    for (int m0 = 0; m0 < inst_.m(); m0++){
        for (int j = 0; j < inst_.n(); j++){
            for (int t=inst_.est(j,m0); t <= inst_.lst(j,m0); t++){
                int var = xIdx_[m0][j][t];
                for (int min = 1; min <= inst_.time(j,m0); min++){
                    // cout << min << " " << inst_.time(j,m0) << " " << j+1 << " " << m0+1 << endl;
                    for (int j1 = 0; j1 < inst_.n(); j1++){
                        //variables_pack.emplace_back();
                        //int pos = variables_pack.size()-1;
                        vector<int> aux_vector;
                        aux_vector.emplace_back(var);
                        //variables_pack[pos].emplace_back(var);
                        // cout << names[var] << " ";
                        for (int j_aux = 0; j_aux < inst_.n(); j_aux++){
                            //cout << "t1: " << min << " " << j_aux+1 << " ";
                            if (j1 == j_aux){
                                for (int t0 = t-inst_.time(j_aux,m0)+1; t0 < t+min; t0++){
                                    if (xIdx_[m0][j_aux][t0] == var) continue;
                                    if (t0 >= inst_.est(j_aux,m0) && t0 <= inst_.lst(j_aux,m0)){
                                        // cout << "[" << j_aux+1 << "," << m0+1 << "," << t0 << "," << t+inst_.minimumTime(m0) << "," << inst_.est(j_aux,m0) << "," << inst_.lst(j_aux,m0) << "] ";
                                        // cout << names[xIdx_[m0][j_aux][t0]] << " ";
                                        // variables_pack[pos].emplace_back(xIdx_[m0][j_aux][t0]);
                                        aux_vector.emplace_back(xIdx_[m0][j_aux][t0]);
                                    }
                                }
                            } else {
                                int qtd_movimentos = min - inst_.time(j_aux,m0);
                                while (qtd_movimentos <= 0){
                                    int t0 = t+qtd_movimentos;
                                    qtd_movimentos++;
                                    if (xIdx_[m0][j_aux][t0] == var) continue;
                                    if (t0 >= inst_.est(j_aux,m0) && t0 <= inst_.lst(j_aux,m0)){
                                        // cout << "(" << j_aux+1 << "," << m0+1 << "," << t0 << "," << t+inst_.minimumTime(m0) << "," << inst_.est(j_aux,m0) << "," << inst_.lst(j_aux,m0) << ") ";
                                        // cout << names[xIdx_[m0][j_aux][t0]] << " ";
                                        // variables_pack[pos].emplace_back(xIdx_[m0][j_aux][t0]);
                                        aux_vector.emplace_back(xIdx_[m0][j_aux][t0]);
                                    }
                                }
                            }
                        }
                        std::sort(aux_vector.begin(), aux_vector.end(), compare);
                        if (aux_vector.size() != 1){
                            unordered_set<vector<int>>::const_iterator got = variables_pack.find(aux_vector);
                            if (got == variables_pack.end()){
                                variables_pack.insert(aux_vector);
                                if (variables_pack.size() % 10000 == 0) cout << j+1 << " " << m0+1 << " " << variables_pack.size() << endl;
                            }
                        }

                    }
                }
            }
        }
    }
    cout << variables_pack.size() << endl;
    ofstream f;
    f.open("cuts.txt");
    for (vector<int> vars : variables_pack){
        for (int v : vars){
            //cout << vars.size() << " " << v << " " ;
            f << names[v] << " ";
        }
        //cout << endl;
        f << endl;
    }
    f.close();
    cout << "arquivo criado" << endl;

    cout << variables_pack.size() << endl;
}


int Fernando::fenchel(int ti, int tf){
    double *x = lp_x(mip);

    unordered_set<vector<S>> solutions;
    vector<S> vars;
    vector<S> vars_valor1;
    for (int t = ti; t < tf; t++){
        for (int j = 0; j < inst_.n(); j++){
            for (int i = 0; i < inst_.m(); i++){
                if (t >= inst_.est(j,i) && t <= inst_.lst(j,i)){
                    // adiciona variáveis possíveis para o corte se o valor relaxado é maior que 0. 
                    // Caso sejam == 1, colocar em outro vetor
                    if (x[xIdx_[i][j][t]] > 0.99999){ 
                        // cout << names[xIdx_[i][j][t]] << " " << x[xIdx_[i][j][t]] << endl;
                        S aux;
                        aux.i = i;
                        aux.j = j;
                        aux.t = t;
                        aux.var = xIdx_[i][j][t];
                        vars_valor1.push_back(aux);
                    } else if (x[xIdx_[i][j][t]] > 1e-05){ 
                        // cout << names[xIdx_[i][j][t]] << " " << x[xIdx_[i][j][t]] << endl;
                        S aux;
                        aux.i = i;
                        aux.j = j;
                        aux.t = t;
                        aux.var = xIdx_[i][j][t];
                        vars.push_back(aux);
                    }
                }
            }
        }
    }

    clock_t begin = clock();
    //int n = vars.size();
    //cout << n << endl;
    vector<S> solution;
    enumeracao_fenchel(vars,0,solutions,solution);

    for (vector<S> v : solutions){
        cout << v.size() << ": ";
        for (S s : v){
            cout << names[s.var] << " ";
        }
        cout << endl;
    }
    getchar();
    if (limite_enumeracao){
        cout << "Limite de enumeração chegado" << endl;
        return 0;
    }
    // cout << "creating set. press enter" << endl;
    // getchar();
    clock_t end = clock();
    cout << "With " << vars.size() + vars_valor1.size() << " variables, found " << solutions.size() << " enumerations in " << (double)(end - begin) / CLOCKS_PER_SEC << " secs" << endl;
    
    // if no enumerations found, end
    if (solutions.size() == 0){
        return 0;
    }
    LinearProgram *fenchel = lp_create();
    lp_set_print_messages(fenchel,0);

    vector< double > lb; // lower bound
    vector< double > ub; // upper bound
    vector< double > obj; // se é objetivo?
    vector< char > integer; // variável inteira?
    vector< string > lambda_names;

    vector< vector< vector< int > > > lambdas = vector<vector<vector<int>>>(inst_.m(),vector<vector<int>>(inst_.n(),vector<int>(inst_.maxTime()+1,-1)));

    unordered_set<S> set_vars;

    for (vector<S> vec : solutions){
        for (S s : vec){
            if (lambdas[s.i][s.j][s.t] == -1){
                lambdas[s.i][s.j][s.t] = lambda_names.size(); 
                lambda_names.push_back("lambda("+to_string(s.j+1)+","+to_string(s.i+1)+","+to_string(s.t)+")"); // nome dessa variável
                lb.push_back(0.0);
                ub.push_back(1.0);
                obj.push_back(-x[s.var]);
                integer.push_back(0);
            }
        }
    }

    lp_add_cols( fenchel, obj, lb, ub, integer, lambda_names );

    int cont = 0;
    for (vector<S> vec : solutions){
        vector< int > idx;
        vector< double > coef;
        for (S s : vec){
            //cout << s.i << " " << s.j << " " << s.t << endl;
            idx.push_back(lambdas[s.i][s.j][s.t]);
            coef.push_back(1);
        }
        lp_add_row(fenchel, idx, coef, "sol("+to_string(cont)+")", 'L', 1.0 );
        cont++;
    }
    
    string filename = inst_.instanceName()+"_fenchel";
    lp_write_lp(fenchel, (filename+".lp").c_str());
    lp_optimize(fenchel);
    
    double *xf = lp_x(fenchel);
    double total = 0;
    for (unsigned int i = 0; i < lambda_names.size(); i++){
        if (xf[i] > 0){
            //cout << lambda_names[i] << " " << xf[i] << endl;

            int mach = 0;
            int job = 0;
            int time = 0;
            sscanf(lambda_names[i].c_str(),"lambda(%d,%d,%d)",&job,&mach,&time);
            
            // cout << lambda_names[i] << "(" << xf[i] <<")" << "*" << names[xIdx_[mach-1][job-1][time]] << "(" << x[xIdx_[mach-1][job-1][time]] <<")" << " + " << endl;
            total += xf[i]*x[xIdx_[mach-1][job-1][time]];
        }
    }
    // cout << total << endl;
    if (total > 1.000001){
        cout << "fenchel violado: " << total << endl;
        vector< int > idx;
        vector< double > coef;
        for (unsigned int i = 0; i < lambda_names.size(); i++){
            if (xf[i] > 0){
                int mach = 0;
                int job = 0;
                int time = 0;
                // cout << lambda_names[i].c_str() << " " << xf[i] << endl;
                sscanf(lambda_names[i].c_str(),"lambda(%d,%d,%d)",&job,&mach,&time);
                idx.push_back(xIdx_[mach-1][job-1][time]);
                coef.push_back(xf[i]);
            }
        }

        for (S s : vars_valor1){
            idx.push_back(xIdx_[s.i][s.j][s.t]);
            coef.push_back(-1.0);
        }
        // for (unsigned int i = 0; i < idx.size(); i++){
        //     cout << coef[i] << "*" << names[idx[i]] << " ";
        // }
        double c = 1 + (double)vars_valor1.size();
        lp_add_row(mip, idx, coef, "fenchel("+to_string(totalFenchel)+")", 'L', c );
        totalFenchel++;
        // lp_write_lp(mip,"teste.lp");
        
    }
    // getchar();
    // delete []x;
    // delete []xf;
    if (total > 1.000001) {
        return 1;
    }
    else{
        return 0;
    }
    //getchar();
}      

// return true if smaller vector is in bigger vector
// return false otherwise
template <typename T> bool Fernando::isSubset(std::vector<T> &A, std::vector<T> &B){
    if (A.size() > B.size()){
        sort(A.begin(), A.end());
        sort(B.begin(), B.end());
        return includes(A.begin(), A.end(), B.begin(), B.end());
    } else {
        sort(A.begin(), A.end());
        sort(B.begin(), B.end());
        return includes(B.begin(), B.end(), A.begin(), A.end());
    }
}

void Fernando::enumeracao_fenchel(const vector<S> &vars, int index, unordered_set<vector<S>> &solutions, vector<S> solution){
    //cout << "solution.size: " << solution.size() << " vars.size: " << vars.size() << " index: " << index << " solutions.size: " << solutions.size() << endl;
    // for (S s : solution){
    //     cout << names[s.var] << " ";
    // }
    // cout << endl;
    // limite de enumerações
    if (solutions.size() > 5000){
        limite_enumeracao = true;
        return;
    }


    // fim da recursão. Adiciona no conjunto de soluções
    if (solution.size() == inst_.m()*inst_.n()){
            // for (S s : solution){
            //     cout << names[s.var] << " ";
            // }
            // cout << endl;
            solutions.insert(solution);
            cout << " solutions: " << solutions.size() << " size solution: " << solution.size() << endl;
        return;
    }
    int qtd = solutions.size();
    // continuação de inserção de variáveis. Caso seja possível inserir, vai pra próxima variável
    for (unsigned int i = index; i < vars.size(); i++){
        if (limite_enumeracao){
            break;
        }
        qtd = solutions.size();
        //cout << solution.size() << " tentando inserir " << names[vars[i].var] << endl;
        //if (insertVar(solution, vars[i])){
        if (canInsert(vars[i])){
            solution.emplace_back(vars[i]);
            enum_time[vars[i].j][vars[i].i] = vars[i].t;
            enumeracao_fenchel(vars,i+1,solutions,solution);
            cout << "qtd_solutions: " << qtd << " solutions: " << solutions.size();
            // no solution inserted
            if (solutions.size() == qtd){
                cout << " size solution: " << solution.size();
                solutions.insert(solution);
            }
            cout << endl;
            solution.pop_back();
            enum_time[vars[i].j][vars[i].i] = -1;
        }
    }
}

// verifica se no conjunto de soluções há alguma que esteja dominando a solução que tentamos inserir.
template <typename T> bool Fernando::dominancia(vector<T> &vec, unordered_set<vector<T>> &set){
    bool exists_solution = false;
    for (vector<T> v : set){
        // caso o tamanho do vetor seja igual, precisamos olhar a dominância apenas baseado na primeira variável.
        // se não for, precisamos olhar independente da ordem
        if (v.size() >= vec.size()){
            exists_solution = isSubset(v,vec);
        } else {
            exists_solution = isSubset(vec,v);
        }
        if (exists_solution){
            return true;
        }
    }
    return false;
}

bool Fernando::insertVar(vector<S> sol, S var){
    int tf_var = var.t + inst_.time(var.j,var.i);
    for (S s : sol){
        if (s.i == var.i || s.j == var.j) {
            int tf_s = s.t + inst_.time(s.j,s.i);

            auto Min = std::max(s.t, var.t);
            auto Max = std::min(tf_s-1, tf_var-1);
            // cout << "Min: " << Min << " Max: " << Max << endl;
            // cout << "insertVar i: " << var.t << " " << tf_var << " " << s.t << " " << tf_s << endl;
            if (Min <= Max) {
                return false;
            }
        }
    }
    return true;
}

int Fernando::executeFenchel(){
    int divisions = 3; // how many time divisions will make
    double interval = inst_.maxTime() / (double)divisions;
    
    interval = teto(interval); // size of time interval
    //cout << interval << endl;
    int qtdCuts = 0;
    for (int i = 0; i < (divisions*2) - 1; i++){
        enum_time = vector< vector< int > >(inst_.n(),vector< int >(inst_.m(),-1));
        limite_enumeracao = false;
        int ti = teto( (i * interval) /2);
        int tf = ti + interval;
        tf = (tf < inst_.maxTime() ? tf : inst_.maxTime());
        cout << "initial: " << ti << " final: " << tf << ". ";
        qtdCuts += fenchel(ti,tf);
        enum_time.clear();
    }
    return qtdCuts;
}

bool Fernando::canInsert(S var){
    int tf_var = var.t + inst_.time(var.j,var.i);

    // check conflict in machine
    for (int i = 0; i < inst_.m(); i++){
        if (var.i == i || enum_time[var.j][i] == -1) continue;
        int tf_s = enum_time[var.j][i] + inst_.time(var.j,i);
        auto Min = std::max(enum_time[var.j][i], var.t);
        auto Max = std::min(tf_s-1, tf_var-1);
        // cout << "Min: " << Min << " Max: " << Max << endl;
        // cout << "insertVar i: " << var.t << " " << tf_var << " " << s.t << " " << tf_s << endl;
        if (Min <= Max) {
            return false;
        }
    }

    // check conflict in job
    for (int j = 0; j < inst_.n(); j++){
        if (var.j == j || enum_time[j][var.i] == -1) continue;
        int tf_s = enum_time[j][var.i] + inst_.time(j,var.i);
        auto Min = std::max(enum_time[j][var.i], var.t);
        auto Max = std::min(tf_s-1, tf_var-1);
        // cout << "Min: " << Min << " Max: " << Max << endl;
        // cout << "insertVar i: " << var.t << " " << tf_var << " " << s.t << " " << tf_s << endl;
        if (Min <= Max) {
            return false;
        }
    }


    int o = inst_.orderMachine(var.j,var.i); // order of the machine;

    if (o != 0) { // not first machine
        // check order conflict
        int m_anterior = inst_.machine(var.j,o-1);
        
        // tarefa anterior foi alocada
        if (enum_time[var.j][m_anterior] != -1){
            int tf_anterior = enum_time[var.j][m_anterior] + inst_.time(var.j,m_anterior);
            auto Min = std::max(enum_time[var.j][m_anterior], var.t);
            auto Max = std::min(tf_anterior-1, tf_var-1);
            // cout << "Min: " << Min << " Max: " << Max << endl;
            // cout << "insertVar i: " << var.t << " " << tf_var << " " << s.t << " " << tf_s << endl;
            if (Min <= Max) {
                return false;
            }        
        }
    }

    if (o != inst_.m()-1){ // not last machine
        int m_posterior = inst_.machine(var.j,o+1);
        
        // tarefa posterior foi alocada
        if (enum_time[var.j][m_posterior] != -1){
            int tf_anterior = enum_time[var.j][m_posterior] + inst_.time(var.j,m_posterior);
            auto Min = std::max(enum_time[var.j][m_posterior], var.t);
            auto Max = std::min(tf_anterior-1, tf_var-1);
            // cout << "Min: " << Min << " Max: " << Max << endl;
            // cout << "insertVar i: " << var.t << " " << tf_var << " " << s.t << " " << tf_s << endl;
            if (Min <= Max) {
                return false;
            }        
        }

    }

    return true;
}