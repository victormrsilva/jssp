#include "Instance.hpp"
#include "lp.hpp"

class Kondili{
public:
    Kondili( const Instance &_inst );

    void optimize();

    virtual ~Kondili();
private:
    const Instance &inst_;

    std::vector< std::vector< std::vector< int > > > xIdx_;

    int cIdx_;
    LinearProgram *mip;
};
