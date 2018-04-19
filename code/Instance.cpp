#include "Instance.hpp"
#include <fstream>
#include <iostream>

using namespace std;

#define str(a) to_string(a)

Instance::Instance( const std::string &fileName, int idx ){
    ifstream ifs;
    ifs.open( fileName.c_str() );

    string line;

    for ( int ii=0 ; ii<=idx ; ++ii )
    {
        getline( ifs, line );
        ifs >> n_ >> m_;

        times_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );
        machines_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );


        // est_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );
        // lst_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );

        cout << "there are " << n() << " jobs and " << m() << " machines." << endl;
        t_ = 0;

        for (int i=0; i < n_; i++){
            int sum_time = 0;
            getline( ifs, line);
            for (int j = 0; j < m_; j++){
                int machine,time;
                ifs >> machine >> time;
                machines_[i][j] = machine;
                times_[i][machine] = time; // para espelhar a formulação matemática do gurobi
                sum_time += time;
                //ifs >> machines_[i][j] >> times_[i][j];                
            }
            if (t_ < sum_time){
                t_ = sum_time;
            }
        }

    //     int seed, machine, lb;
    //     string tmp;
    //     ifs >> seed >> machine >> lb >> ub_ >> tmp;

    //     for ( int i=0 ; (i<n_) ; ++i )
    //         for ( int j=0 ; (j<m_) ; ++j )
    //             ifs >> times_[i][j];

    //     ifs >> tmp;


    //     for ( int i=0 ; (i<n_) ; ++i )
    //     {
    //         for ( int j=0 ; (j<m_) ; ++j )
    //         {
    //             ifs >> machines_[i][j];
    //             --machines_[i][j];
    //         }
    //     }
    //     getline( ifs, line );
    // }

    // // computing est and lst
    // for ( int j=0 ; (j<n_) ; ++j )
    // {
    //     int t=0;
    //     for ( int i=0 ; (i<m_) ; ++i )
    //     {
    //         int ii = machines_[j][i];
    //         est_[j][ii] = t;
    //         t += times_[j][ii];
    //     }
    // }

    // for ( int j=0 ; (j<n_) ; ++j )
    // {
    //     int t=ub_;
    //     for ( int i=m_-1 ; (i>=0) ; --i )
    //     {
    //         int ii = machines_[j][i];
    //         t -= times_[j][ii];
    //         lst_[j][ii] = t;
    //     }
    }




    ifs.close();
}

void Instance::saveCmpl( const string fname ) const {
    ofstream out( fname );

    out << "%Jobs set <0.."+str(n()-1)+">" << endl;;
    out << "%Machines set <0.."+str(n()-1)+">" << endl << endl;

    out << "%order[Jobs,Machines] < ";

    for ( int j=0 ; (j<n()) ; ++j )
    {
        for ( int i=0 ; (i<n()) ; ++i )
            out << machine( j, i ) << " ";
        if (j<n()-1) 
            out << endl << "    ";
        else
            out << " >" << endl;
    }
    out << endl;

    out << "%time[Jobs,Machines] < ";

    for ( int j=0 ; (j<n()) ; ++j )
    {
        for ( int i=0 ; (i<n()) ; ++i )
            out << time( j, i ) << " ";
        if (j<n()-1) 
            out << endl << "    ";
        else
            out << " >" << endl;
    }
    out << endl;

    out.close();
}

