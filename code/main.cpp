#include <string>
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <numeric>
#include <time.h>
#include "src/Instance.hpp"
#include "src/Compact.hpp"
#include "src/Flow_testes.hpp"
#include "src/Flow.hpp"
#include "src/Kondili.hpp"
#include "src/Fernando.hpp"
#include "src/Gera.hpp"


using namespace std;
/**
 * @brief initialization: jssp instanceFile time. instanceFile is the name of the instance. time is the time of reference. 0 means no time and thus it will be calculated
 * 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */

constexpr unsigned int str2int(const char* str, int h = 0){
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

int RandomNumber () { return (std::rand()%10)+1; }

void geraInstancias( string filename, int tamanho){
    vector<int> machines(tamanho);
    vector<int> time(tamanho);
    iota(machines.begin(),machines.end(),0);
    ofstream out(filename);
    
    out << tamanho << " " << tamanho << endl;
    for (int j = 0; j < tamanho; j++){
        random_shuffle ( machines.begin(), machines.end() );
        generate(time.begin(), time.end(), RandomNumber);
        for (int k = 0; k < tamanho; k++){
            out << machines[k] << " " << time[k] << " " ;
        }
        out << endl;
    }
    out.close();
}

int main( int argc, char **argv )
{
    srand(time(NULL));
    if (argc < 5)
    {
        cerr << "jssp instanceFile timeLimit execute formulation" << endl;
        cerr << "timeLimit: a integer with limited time of execution. -1 will be maximum time needed" << endl;
        cerr << "execute: 1 for execute the lp. 0 only generates lp" << endl;
        cerr << "formulation: " << endl;
        cerr << "F for Flow " << endl;
        cerr << "C for Compact (BigM) " << endl;
        cerr << "K for Kondilli " << endl;
        cerr << "Fe for Fernando " << endl;
        cerr << "T for tests " << endl;
        exit(EXIT_FAILURE);
    }



    

    string option = string(argv[4]);

    if (option == "G"){
        double maximo = 0;
        int max = 0;
        for (int i = 0; i < 50; i++){
            cout << i << ":" << endl;
            string filename = "vi"+to_string(i);
            geraInstancias(filename,3);
            Instance inst(filename,90,1);
            Gera mip( inst );
            double valor = mip.execute();
            if (valor > maximo){
                maximo = valor;
                max = i;
            }
            //getchar();
        }
        cout << "Maximo: " << max << " " << maximo << endl;
        return EXIT_SUCCESS;
    }
    Instance inst( string(argv[1]), atoi(argv[2]), atoi(argv[3]) );

    inst.saveCmpl("jssp.cdat");
    cout << inst.m() << " " << inst.n() << endl;
    if (option == "F"){
        Flow mip( inst );
    }
    if (option == "C"){
        Compact mip( inst );
    }

    if (option == "K"){
        Kondili mip( inst );
    }

    if (option == "Fe"){
        cout << "Fernando formulation selected" << endl;
        Fernando mip( inst );
    }

    if (option == "T"){
        cout << "Test formulation selected" << endl;
        Flow_testes mip( inst );
        mip.teste_elimina_variavel();
    }
    

    return EXIT_SUCCESS;
}
