#include "Instance.hpp"
#include "lp.hpp"

class Fernando{
public:
    Fernando( const Instance &_inst );

    void optimize();

    virtual ~Fernando();
private:
    const Instance &inst_;

    std::vector< std::vector< std::vector< int > > > xIdx_;
    std::vector< std::vector< std::vector< int > > > eIdx_;
    std::vector< std::vector< int > > fIdx_;

    int cIdx_;
    LinearProgram *mip;
};
