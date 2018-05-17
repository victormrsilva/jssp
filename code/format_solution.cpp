#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

vector<string> explode(string const & s, char delim)
{
    vector<string> result;
    istringstream iss(s);

    for (string token; getline(iss, token, delim); )
    {
        result.push_back(move(token));
    }

    return result;
}

vector< vector< int > > times;

vector< vector< int > > machines;
vector< vector< int > > ordem;

// earliest starting time of job on machine
vector< vector< int > > est_;

// latest starting job of job on machine
vector< vector< int > > lst_;

int n,m;

int main(int argc, char** argv){
  if (argc < 3){
    cout << "format_solution instance_file solution_file new_solution_file";
  }

  ifstream ifs;
  ofstream ofs;
  ifs.open( argv[1] );

  string line;

  //getline( ifs, line );
  ifs >> n >> m;

  times = vector< vector< int > >( n, vector<int>( m, 0 ) );
  machines = vector< vector< int > >( n, vector<int>( m, 0 ) );
  ordem = vector< vector< int > >( n, vector<int>( m, 0 ) );


  est_ = vector< vector< int > >( n, vector<int>( m, 0 ) );
  lst_ = vector< vector< int > >( n, vector<int>( m, 0 ) );

  cout << "there are " << n << " jobs and " << m << " machines." << endl;

  for (int i=0; i < n; i++){
      getline( ifs, line);
      for (int j = 0; j < m; j++){
        int machine,time;
        ifs >> machine >> time;
        machines[i][j] = machine;
        ordem[i][machine] = j;
        times[i][machine] = time; // para espelhar a formulação matemática do gurobi
        //cout << i << " " << j <<" " << machines[i][j] << " " <<  times[i][machine] << " " << ordem[i][machine] << endl;
      }
  }
  


  ifs.close();

  ifs.open( argv[2] );

  ofs.open( argv[3] );

  while (getline(ifs, line)){
    if (line.find('x') != string::npos){ // encontrou x
      // quebrar solução. primeiro pegar o valor
      vector<string> exploded = explode(line,' ');
      int tempo = stoi(exploded[1]);

      exploded = explode(exploded[0],'(');
      exploded = explode(exploded[1],')');
      exploded = explode(exploded[0],',');

      int job = stoi(exploded[0]);
      int machine = stoi(exploded[1]);
      int proxima = ordem[job-1][machine-1] + 1;
      //cout << job << " " << machine << " " << tempo << endl;
      //cout << " x(" << job << "," << machine << "," << tempo << "," << (proxima == m ? proxima+1 : machines[job-1][proxima]+1 ) << "," << tempo+times[job-1][machine-1] << ")" << endl;
      string var  = " x(" + to_string(job) + "," + to_string(machine) + "," + to_string(tempo) + "," + (proxima == m ? "f" : to_string(machines[job-1][proxima]+1) ) + "," + to_string(tempo+times[job-1][machine-1]) + ")";
      //replace( var, ","+to_string(m+1)+",", ",f," );
      
      ofs << var << " " << 1 << endl;
      cout << var << " " << 1 << endl;
    }
  }
    ifs.close();
    ofs.close();
}