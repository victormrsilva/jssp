#include "Instance.hpp"
#include <fstream>
#include <iostream>

using namespace std;

#define str(a) to_string(a)

Instance::Instance( const std::string &fileName, int time , int execute){
    ifstream ifs;
    ifs.open( fileName.c_str() );

    instance_name = fileName.c_str();

    execute_ = execute;

    string line;

    //getline( ifs, line );
    n_ = m_ = 0;
    ifs >> n_ >> m_;

    times_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );
    minimum_time_ =  vector< int >( m_, __INT_MAX__);
    distances_ = vector< vector< vector < int > > >( n_, vector< vector <int> >( m_, vector<int>(m_,0)) );
    machines_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );
    order_ = vector< vector< int > >( n_, vector<int>( m_, -1 ) );


    est_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );
    lst_ = vector< vector< int > >( n_, vector<int>( m_, 0 ) );

    cout << "there are " << n() << " jobs and " << m() << " machines." << endl;
    bool given_time = false;
    if (time == 0){
        t_ = 0; // worse time of completion (all jobs happening simultaniously)
        given_time = false;
    } else {
        t_ = time;
        given_time = true;
    }
        
    

    for (int i=0; i < n_; i++){
        getline( ifs, line);
        for (int j = 0; j < m_; j++){
            int machine,time;
            ifs >> machine >> time;
            machines_[i][j] = machine;
            order_[i][machine] = j;
            times_[i][machine] = time; // para espelhar a formulação matemática do gurobi
            if (time < minimum_time_[machine]){
                minimum_time_[machine] = time;
            }
            if (!given_time) t_ += time;
        }
    }
    

    // computing distance
    for ( int j=0 ; (j<n_) ; ++j ){
        for ( int i=0 ; (i<m_) ; ++i ){
            int machine_a = machines_[j][i];
            for ( int k=i+1 ; (k<m_) ; ++k ){
                int machine_b = machines_[j][k];
                int dur = times_[j][machines_[j][i]];
                for (int aux = i+1; aux < k; aux++){
                    dur += times_[j][machines_[j][aux]];
                }
                distances_[j][machine_a][machine_b] += dur ;
            }
        }
    }    


    // computing est and lst
    for ( int j=0 ; (j<n_) ; ++j ){
        int t=0;
        for ( int i=0 ; (i<m_) ; ++i ){
            int mach = machines_[j][i]; // machine in order i for job j
            est_[j][mach] = t;
            t += times_[j][mach];
        }
    }

    for ( int j=0 ; (j<n_) ; ++j ){
        int t=t_;
        for ( int i=m_-1 ; (i>=0) ; --i ){ // start from last
            int mach = machines_[j][i];
            t -= times_[j][mach];
            lst_[j][mach] = t;
        }
    }

    ifs.close();
}

int Instance::execute() const {
    return execute_;
}

void Instance::saveCmpl( const string fname ) const {
    ofstream out( fname );

    out << "%Jobs set <0.."+str(n()-1)+">" << endl;;
    out << "%Machines set <0.."+str(m()-1)+">" << endl << endl;

    out << "%order[Jobs,Machines] < " << endl;

    for ( int j=0 ; (j<n()) ; ++j )
    {
        for ( int i=0 ; (i<m()) ; ++i )
            out << machine( j, i )+1 << " ";
        if (j<n()-1) 
            out << endl;
        else
            out << " >" << endl;
    }
    out << endl;

    out << "%time[Jobs,Machines] < " << endl;

    for ( int j=0 ; (j<n()) ; ++j )
    {
        for ( int i=0 ; (i<m()) ; ++i )
            out << time( j, i ) << " ";
        if (j<n()-1) 
            out << endl;
        else
            out << " >" << endl;
    }
    out << endl;
    out << "MaxTime: " << maxTime() << endl;
    out << "Job n, Machine m: est lst" << endl;
    for ( int j=0 ; (j<n()) ; ++j )    {
        for ( int i=0 ; (i<m()) ; ++i ){
            out << "Job " << j+1 << ", Machine " << i+1 <<": " << est(j,i) << " " << lst(j,i) << endl;
        }
    }
    out << endl;
    out << "Job n, Machine i -> Machine j: distance" << endl;
    for ( int j=0 ; (j<n()) ; ++j )    {
        for ( int i=0 ; (i<m()) ; ++i ){
            for ( int k=i+1 ; (k<m()) ; ++k ){
                out << "Job " << j+1 << ", Machine " << machines_[j][i]+1 <<" -> Machine " << machines_[j][k]+1 << ": " << distance(j,machines_[j][i], machines_[j][k])  << endl;
            }
        }
    }
    out << endl;
    out << "Machine i: minimum_time " << endl;
    for ( int j=0 ; (j<n()) ; ++j )    {
        out << "Machine " << j+1 << ": " << minimum_time_[j] << endl;
    }

    out.close();
}

string Instance::instanceName() const{
    return instance_name;
}

