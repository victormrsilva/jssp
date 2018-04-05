#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include <string>
#include <vector>

/**
 * @brief class that generates the instance
 * 
 */

class Instance
{
public:

    /**
     * @brief Construct a new Instance object
     * 
     * The instance must be in form (machine start with 0)
     * number_jobs number_machines
     * job1_machine job1_time_machine job1_machine job1_time_machine ... job1_machine job1_time_machine (m times)
     * job2_machine job2_time_machine job2_machine job2_time_machine ... job2_machine job2_time_machine (m times)
     * ...
     * jobn_machine jobn_time_machine jobn_machine jobn_time_machine ... jobn_machine jobn_time_machine (m times)
     * 
     * As described in http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/jobshop1.txt
     * 
     * @param fileName instance filename. 
     * @param idx number of lines in instance
     * 
     */
    Instance( const std::string &fileName, int idx = 0 );

    /**
     * @brief return number of jobs
     * 
     * @return int 
     */
    int n() const { return n_; } // number of jobs

    /**
     * @brief return number of machines
     * 
     * @return int 
     */
    int m() const { return m_;} // number of machines

    /**
     * @brief processing time of job j on machine i
     * 
     * @param j - job. Jobs starts with -
     * @param i - machine
     * @return int 
     */
    int time( int j, int i ) const { return times_[j][i]; }

    /**
     * @brief order of machines for the job, which is the i-the machine where job j must be processed
     * 
     * @param j - job
     * @param i - machine
     * @return int 
     */
    int machine( int j, int i ) const { return machines_[j][i]; }

    // known upper bound
    int ub() const { return ub_; }

    // earliest starting time of job on machine
    int est( int j, int i ) const { return est_[j][i]; }

    // latest starting job of job on machine
    int lst( int j, int i ) const { return lst_[j][i]; }

    // save the cmpl file
    void saveCmpl( const std::string fname ) const;
private:
    int n_;
    int m_;

    int ub_;

    std::vector< std::vector< int > > times_;

    std::vector< std::vector< int > > machines_;

    // earliest starting time of job on machine
    std::vector< std::vector< int > > est_;

    // latest starting job of job on machine
    std::vector< std::vector< int > > lst_;
};

#endif

