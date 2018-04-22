#include <string>
#include <iostream>
#include "Instance.hpp"
//#include "Compact.hpp"
#include "Flow.hpp"

using namespace std;
/**
 * @brief initialization: jssp instanceFile time. instanceFile is the name of the instance. time is the time of reference. 0 means no time and thus it will be calculated
 * 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main( int argc, char **argv )
{
    if (argc<3)
    {
        cerr << "jssp instanceFile instanceNum" << endl;
        exit(EXIT_FAILURE);
    }

    Instance inst( string(argv[1]), atoi(argv[2]) );

    inst.saveCmpl("jssp.cdat");

    cout << inst.m() << " " << inst.n() << endl;

    Flow mip( inst );

    return EXIT_SUCCESS;
}
