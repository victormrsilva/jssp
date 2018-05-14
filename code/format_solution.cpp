#include <iostream>
#include <string>
#include <vector>

using namespace std;

vector< vector< int > > times_;

vector< vector< int > > machines_;

// earliest starting time of job on machine
vector< vector< int > > est_;

// latest starting job of job on machine
vector< vector< int > > lst_;

int n,m;

int main(int argc, char* argv){{
  if (argc < 3){
    cout << "format_solution instance_file solution_file new_solution_file"
  }

  ifstream ifs;
  ifs.open( argv[1] );

  string line;

  //getline( ifs, line );
  ifs >> n >> m;

  times_ = vector< vector< int > >( n, vector<int>( m, 0 ) );
  machines_ = vector< vector< int > >( n, vector<int>( m, 0 ) );


  est_ = vector< vector< int > >( n, vector<int>( m, 0 ) );
  lst_ = vector< vector< int > >( n, vector<int>( m, 0 ) );

  cout << "there are " << n << " jobs and " << m << " machines." << endl;

  for (int i=0; i < n; i++){
      getline( ifs, line);
      for (int j = 0; j < m; j++){
        int machine,time;
        ifs >> machine >> time;
        machines_[i][j] = machine;
        times_[i][machine] = time; // para espelhar a formulação matemática do gurobi
      }
  }
  

  // computing est and lst
  for ( int j=0 ; (j<n) ; ++j ){
      int t=0;
      for ( int i=0 ; (i<m) ; ++i ){
        int mach = machines_[j][i]; // machine in order i for job j
        est_[j][mach] = t;
        t += times_[j][mach];
      }
  }

  for ( int j=0 ; (j<n) ; ++j ){
      int t=t_;
      for ( int i=m-1 ; (i>=0) ; --i ){ // start from last
        int mach = machines_[j][i];
        t -= times_[j][mach];
        lst_[j][mach] = t;
      }
  }

  ifs.close();


}