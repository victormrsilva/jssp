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

    std::vector< std::vector< std::vector< std::vector< std::vector< int > > > > > xIdx_; // índice dos arcos (y,m0,t0,mf,tf)

    //std::vector< std::vector< std::vector< int > > > yIdx_;

    int cIdx_;
    LinearProgram *mip;

    // std::vector < std::string > x11_color_;

    // void create_x11();
    void createCompleteGraphDot();
};
