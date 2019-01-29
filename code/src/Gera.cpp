

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





Gera::Gera( const Instance &_inst ) : inst_(_inst),
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

    //unordered_set<std::vector<int>, VectorHash> variables_pack;

    
    
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
                        bool vetor_igual = false;
                        std::sort(aux_vector.begin(), aux_vector.end(), [](int a, int b) {return a < b;});
                        if (aux_vector.size() != 1){
                            unordered_set<vector<int>>::const_iterator got = variables_pack.find(aux_vector);
                            if (got == variables_pack.end()){
                                variables_pack.insert(aux_vector);
                                if (variables_pack.size() % 10000 == 0) cout << j+1 << " " << m0+1 << " " << variables_pack.size() << endl;
                            }
                        }
                        
                        // for (unsigned v=0; v < variables_pack.size(); v++){
                        //     if (variables_pack[v] == aux_vector){
                        //         vetor_igual = true;
                        //         v=variables_pack.size();
                        //     }
                        // }
                        // if (!vetor_igual && aux_vector.size() != 1){
                        //     variables_pack.emplace_back(aux_vector);
                        // }
                        // cout << endl;
                    }
                }
            }
        }
    }
    
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
    //getchar();
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
    // if (inst_.execute()){
    //     lp_optimize( mip );
    //     lp_write_sol(mip, "jssp_Fernando.sol");
    // }

}

int Gera::manual_cuts(){
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
    string filename = inst_.instanceName()+"_cortes";
    lp_write_lp(mip, (filename+".lp").c_str());
    //getchar();
    return qtd;
}

double Gera::execute(){
    lp_set_print_messages(mip,0);
    lp_write_lp(mip,"teste.lp");        
    lp_optimize(mip);
    double bnd_integer = lp_obj_value(mip);
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
    
    lp_optimize_as_continuous(mip);
    lp_write_lp(mip,"teste_cb.lp");
    double bnd_continuous = lp_obj_value(mip);
    cout << "Continuo = " << bnd_continuous << endl;
    lp_write_lp(mip,"teste_cb.lp");

    lp_as_integer(mip);

    //Callback cb = Callback(mip,inst_,xIdx_,process);
    
    //CGraph *cgraph = build_cgraph_lp(mip);
    //cgraph_print_summary(cgraph, "test_cgraph");


    cout << "Inteiro = " << bnd_integer << endl;
    cout << "Diferença = " << bnd_integer - bnd_continuous << endl;
    //getchar();
    //lp_write_lp(mip,"teste_cb.lp");
    lp_write_sol(mip,"teste_cb.sol");
    return bnd_integer - bnd_continuous;
}

Gera::~Gera()
{
    lp_free( &mip );
}



double Gera::lifting(double c, int *idxs, double *coefs){
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
        //cout << "C: " << c << " soma: " << sol << " soma - C: " << violado << endl;
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

void Gera::lifting_binario(int *idxs, double *coefs){
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
        
        bnd = lifting(teto(c),&idxs[0],&coefs[0]);
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

void Gera::lifting_linear(int *idxs, double *coefs){

    clock_t begin = clock();
    int iteracoes = 0;
    lp_optimize_as_continuous(mip);
    double bnd = lp_obj_value(mip);
    double bnd_anterior = 0;
    double c = 0;
    ofstream saida(inst_.instanceName()+"_solution_linear_machine.csv");
    saida << "iteracao;bnd;c;cortes;tempo" << endl;
    saida << "0;" << bnd << ";0;0;"<<lp_solution_time(mip)<<endl;
    int cuts = 1;
    // getchar();
    // while (cuts != 0) {
    //     cout << "do" << endl;
    //     lp_write_lp(mip,"teste.lp");
    //     lp_optimize_as_continuous(mip);
    //     cout << "do1" << endl;
    //     cuts = oddHoles();
    //     cout << "oddholes: " << cuts << endl;
    //     getchar();
    //     // if (cuts == 0){
    //     //     cuts = cliques(idxs,coefs);
    //     //     cout << "cortes cliques: " << cuts << endl;
    //     //     getchar();
    //     // }
    // };
    while (fabs(bnd - teto(bnd_anterior)) > 1e-06) {
        iteracoes++;
        int cortes = 0;
        if (clique)
        {
            cortes = manual_cuts();
            //cortes = oddHoles();
            //cout << "oddHoles: " << cortes << endl;
            //cortes = cortes + cliques(idxs,coefs);
        }
        
        bnd_anterior = bnd;
        c = teto(bnd);
        bnd = lifting(c,idxs,coefs);
        
        //lp_write_lp(mip, "teste_cb.lp");
        //lp_write_sol(mip, "solution.sol");
        saida << iteracoes << ";" << bnd << ";" << c << ";" << cortes << ";" << lp_solution_time(mip)<<endl;
        //cout <<  " c: " << c << " bnd_anterior: " << bnd_anterior << " bnd: " << bnd << endl;
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