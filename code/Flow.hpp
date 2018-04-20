#include "Instance.hpp"
#include "lp.hpp"
class Flow
{
public:
    Flow( const Instance &_inst );

    void optimize();

    virtual ~Flow();
private:
    const Instance &inst_;

    std::vector< std::vector< std::map< std::map< std::map< int > > > > > xIdx_; // índice dos arcos (y,m0,t0,mf,tf)

    //std::vector< std::vector< std::vector< int > > > yIdx_;

    int cIdx_;
    LinearProgram *mip;

    int getXidx(int j, int m0, int t0, int mf, int tf) const;

    // std::vector < std::string > x11_color_;

    // void create_x11();
    void createCompleteGraphDot();

      // set of entering flows
    vector<vector<vector<int>>> enter_flow;
    // set of exiting flows
    vector<vector<vector<int>>> exit_flow;
    // proccessing jobs flows
    vector<vector<vector<int>>> proccess;

};
