#include "Instance.hpp"
#include "lp.hpp"

class Compact
{
public:
    Compact( const Instance &_inst );

    void optimize();

    virtual ~Compact();
private:
    const Instance &inst_;

    std::vector< std::vector< int > > xIdx_;

    std::vector< std::vector< std::vector< int > > > yIdx_;

    int cIdx_;
    LinearProgram *mip;
};
