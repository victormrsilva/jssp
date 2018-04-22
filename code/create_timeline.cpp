#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

using namespace std;

int main( int argc, char **argv )
{
    if (argc<3)
    {
        cerr << "create_timeline machine maxTime" << endl;
        exit(EXIT_FAILURE);
    }

    
    int m = stoi(argv[1]);
    int time = stoi(argv[2]);
    ifstream file("solution.txt");
    string line;
    int num,j,m0,t0,mf,tf,one,zero;

    vector<vector<string>> timeline;

    timeline = vector<vector<string>>(time+1, vector<string>(m+2));


    while( getline( file, line ) )
    {
        
        sscanf(line.c_str(), "%d x(%d,%d,%d,%d,%d) %d %d", &num, &j, &m0, &t0, &mf, &tf,&one,&zero);
        cout << j << " " << m0 << " " << t0 << " " << mf << " " << tf << endl;
        if (m0 == mf){
            timeline[t0][m0]= timeline[t0][m0]+","+to_string(j)+"w";
        } else {
            for (int i = t0; i < tf; i++ ){
                timeline[i][m0] = timeline[i][m0]+","+to_string(j)+"p";
            }
        }
    }

    ofstream f;
    f.open ("timeline.txt");
    f << "-";
    for (int i = 1; i <= m; i++){
        f << "\t " << i;
    }
    f << endl;
    for (int t = 1; t < time+1; t++){
        f << t;
        for (int m0 = 1; m0 <= m; m0++){
            f << "\t" << timeline[t][m0];
        }
        f << endl;
    }
    f.close();

}