#include <string>
#include <iostream>
#include "src/Instance.hpp"
#include "src/Compact.hpp"
#include "src/Flow.hpp"
#include "src/Kondili.hpp"


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

int main( int argc, char **argv )
{
    if (argc < 4)
    {
        cerr << "jssp instanceFile instanceNum formulation" << endl;
        cerr << "formulation: " << endl;
        cerr << "F for Flow " << endl;
        cerr << "C for Compact (BigM) " << endl;
        exit(EXIT_FAILURE);
    }

    Instance inst( string(argv[1]), atoi(argv[2]) );

    inst.saveCmpl("jssp.cdat");

    cout << inst.m() << " " << inst.n() << endl;

    string option = string(argv[3]);

    if (option == "F"){
        Flow mip( inst );
    }
    if (option == "C"){
        Compact mip( inst );
    }

    if (option == "K"){
        Kondili mip( inst );
    }
    

    return EXIT_SUCCESS;
}
