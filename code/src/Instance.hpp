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
    Instance( const std::string &fileName, int idx = 0, int execute = 0 );

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

    int orderMachine(int job, int machine) { return order_[job][machine]; }
    int distance( int j, int i, int k ) const { return distances_[j][i][k]; }

    // known upper bound
    int ub() const { return ub_; }

    // max time
    int maxTime() const { return t_; };

    // earliest starting time of job on machine i
    int est( int j, int i ) const { return est_[j][i]; }

    // latest starting job of job on machine i
    int lst( int j, int i ) const { return lst_[j][i]; }

    // minimum processing time of machine i
    int minimumTime( int i ) const { return minimum_time_[i]; }

    void setEst(int j, int i, int value) { est_[j][i] = value;}
    void setLst(int j, int i, int value) { lst_[j][i] = value;}

    // save the cmpl file
    void saveCmpl( const std::string fname ) const;

    // return instance name
    std::string instanceName() const;

    // return if will be executed
    int execute() const; 

    // return a copy of est
    std::vector< std::vector <int> > copyEst() { return est_; }
    // return a copy of lst
    std::vector< std::vector <int> > copyLst() { return est_; }
private:
    int n_;
    int m_;

    int t_; // max time

    int ub_;

    int execute_;

    std::string instance_name;

    std::vector< std::vector< int > > times_;
    std::vector< std::vector< std::vector <int > > > distances_;
    std::vector< std::vector< int > > order_;

    std::vector< std::vector< int > > machines_;

    // earliest starting time of job on machine
    std::vector< std::vector< int > > est_;

    // latest starting job of job on machine
    std::vector< std::vector< int > > lst_;

    // minimum processing time of machines
    std::vector<int> minimum_time_;
};

#endif

