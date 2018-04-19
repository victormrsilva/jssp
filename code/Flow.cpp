#include "Flow.hpp"

#include <vector>
#include <string>
#include <cfloat>

using namespace std;

Flow::Flow( const Instance &_inst ) :
    inst_(_inst)
{
    mip = lp_create();

    vector< string > names;
    vector< double > lb;
    vector< double > ub;
    vector< double > obj;
    vector< char > integer;

    xIdx_ = vector< vector< int > >( inst_.n(), vector<int>( inst_.m(), -1 ) );

    yIdx_ = vector< vector< vector<int> > >( inst_.n(), vector< vector< int > >( inst_.n(),
                vector<int>( inst_.m(), -1 ) ) );

    // creating x vars
    for ( int j=0 ; (j<inst_.n()) ; ++j )
    {
        for ( int i=0 ; (i<inst_.m()) ; ++i )
        {
            xIdx_[j][i] = names.size();
            names.push_back( "x("+to_string(j+1)+","+to_string(i+1)+")" );
            lb.push_back( 0.0 );
            ub.push_back( inst_.lst( j, i ) );
            obj.push_back( 0.0 );
            integer.push_back( 1 );
        }
    }

    // y vars
    for ( int j1=0 ; (j1<inst_.n()) ; ++j1 )
    {
        for ( int j2=j1+1 ; (j2<inst_.n()) ; ++j2 )
        {
            for ( int i=0 ; (i<inst_.m()) ; ++i )
            {
                yIdx_[j1][j2][i] = names.size();
                names.push_back( "y("+to_string(j1+1)+","+to_string(j2+1)+","+to_string(i+1)+")" );
                lb.push_back( 0.0 );
                ub.push_back( 1.0 );
                obj.push_back( 0.0 );
                integer.push_back( 1 );
             }
        }
    }

    // c var
    cIdx_ = names.size();
    names.push_back("C");
    lb.push_back( 0.0 );
    ub.push_back( DBL_MAX );
    obj.push_back( 1.0 );
    integer.push_back( 1 );

    lp_add_cols( mip, obj, lb, ub, integer, names );

    // constraint for job on each one of its machines
    for ( int j=0 ; (j<inst_.n()) ; ++j )
    {
        for ( int i=1 ; (i<inst_.m()) ; ++i )
        {
            vector< int > idx;
            vector< double > coef;

            idx.push_back( xIdx_[j][inst_.machine(j,i)] );
            coef.push_back( 1.0 );
            idx.push_back( xIdx_[j][inst_.machine(j,i-1)] );
            coef.push_back( -1.0 );

            lp_add_row( mip, idx, coef, "prec("+to_string(j+1)+","+to_string(i+1)+")", 'G', inst_.time( j, inst_.machine(j,i-1)) );
        }
    }

    // linking c and x
    for ( int j=0 ; (j<inst_.n()) ; ++j )
    {
        vector< int > idx;
        vector< double > coef;

        idx.push_back( cIdx_ );
        coef.push_back( 1.0 );

        idx.push_back( xIdx_[j][inst_.machine(j, inst_.m()-1)] );
        coef.push_back( -1.0 );

        lp_add_row( mip, idx, coef, "lnkCX("+to_string(j+1)+")", 'G', inst_.time( j, inst_.machine(j, inst_.m()-1)) );
    }

    for ( int j1=0 ; (j1<inst_.n()) ; ++j1 )
    {
        for ( int j2=0 ; (j2<inst_.n()) ; ++j2 )
        {
            if (j1==j2)
                continue;
            
            for ( int i=0 ; (i<inst_.m()) ; ++i )
            {
                vector< int > idx; vector< double > coef;
                idx.push_back( xIdx_[j1][i] );
                coef.push_back( 1.0 );
                idx.push_back( xIdx_[j2][i] );
                coef.push_back( -1.0 );
    
                double rhs = inst_.time( j2, i );

                if (j1<j2)
                {
                    rhs -= 99999;
                    idx.push_back( yIdx_[j1][j2][i] );
                    coef.push_back( -99999 );
                }
                else
                {
                    idx.push_back( yIdx_[j2][j1][i] );
                    coef.push_back( 99999 );
                }

                lp_add_row( mip, idx, coef, "lnkXY("+to_string(j1+1)+","+to_string(j2+1)+","+to_string(i+1)+")", 'G', rhs );
            }
        }
    }
        
    lp_optimize( mip );
    lp_write_lp( mip, "jssp" );
}

Flow::~Flow()
{
    lp_free( &mip );
}

