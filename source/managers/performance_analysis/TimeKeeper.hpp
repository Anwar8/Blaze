/**
 * @file TimeKeeper.hpp
 * @brief defines the TimeKeeper class which keeps track of and controls a set of ExecutionTimers.
 */

#ifndef TIME_KEEPER
#define TIME_KEEPER
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include "ExecutionTimer.hpp"

/**
 * @brief Keeps track of and controls a set of \ref ExecutionTimers. Use for profiling.
 * 
 */
class TimeKeeper 
{
    protected:
        std::map<std::string, ExecutionTimer> timers_map; /**< a map of \ref ExecutionTimers that the TimeKeeper is in charge of.*/
        int rank = 0;
        int num_ranks = 1;
        std::map<int, std::vector<double>> rank_durations_map; /**< a map of ranks pointing to vectors of times kept on all processes.*/
    public:
        /**
         * @brief sets the rank and num_ranks of the time keeper.
         * 
         * @param rank rank of the calling MPI process.
         * @param num_ranks number of available MPI processes.
         */
        void initialise_parallel_keeper(int rank, int num_ranks)
        {
            this->rank = rank;
            this->num_ranks = num_ranks;
        }
        /**
         * @brief populates the \ref timers_map with \ref ExecutionTimer objects that this function will create.
         * 
         * @param timers_names a vector containing the names of each \ref ExecutionTimer object.
         */
        void add_timers(std::vector<std::string> timers_names)
        {
            for (std::string name : timers_names)
            {
                timers_map.insert({name, ExecutionTimer()});
            }
        }

        /**
         * @brief starts a particular timer by calling \ref ExecutionTimer::start.
         * 
         * @param timer_name the timer in the \ref timers_map which will be started.
         */
        void start_timer(std::string timer_name)
        {
            timers_map[timer_name].start();
        }
        
        /**
         * @brief stops a particular timer by calling \ref ExecutionTimer::stop.
         * 
         * @param timer_name the timer in the \ref timers_map which will be stopped.
         */
        void stop_timer(std::string timer_name)
        {
            timers_map[timer_name].stop();
        }
        
        /**
         * @brief resets a particular timer by calling \ref ExecutionTimer::reset.
         * 
         * @param timer_name the timer in the \ref timers_map which will be reset.
         */
        void reset_timer(std::string timer_name)
        {
            timers_map[timer_name].reset();
        }

        /**
         * @brief returns the duration value of a timer by calling its \ref ExecutionTimer::get_duration.
         * 
         * @param timer_name the timer in the \ref timers_map which duration will be retrieved.
         * @return duration of `timer_name`.
         */
        double get_timer_duration(std::string timer_name)
        {
            return timers_map[timer_name].get_duration();
        }

        /**
         * @brief prints a table with durations recorded by all requested timers as well as a percentage calculated with respect to a reference_timer.
         * @details the percentage is calculated using: \f$ P_i = \frac{\Delta t_i \times 100}{\Delta t_{reference}} \f$
         * @param timers_names the names of the timers which durations are to be printed.
         * @param reference_timer the timer which duration will be used to calculate percentages.
         */
        void read_timers(std::vector<std::string> timers_names, std::string reference_timer = "")
        {
            
            std::cout << std::endl <<  "Reading some ExecutionTimers on rank" << rank << ":" << std::endl;
            std::cout << std::left << std::setw(32) << "Timer Name" << ": "
                  << std::setw(12) << "Duration (s)";
            if (reference_timer != "")
            {
                std::cout << std::right << std::setw(12) << "\%"<< std::endl;
                std::cout << std::string(57, '-') << std::endl;

            }
            else
            {
                std::cout << std::endl;
                std::cout << std::string(46, '-') << std::endl;
            }
            
            for (std::string name : timers_names)
            {
                std::cout << std::left << std::setw(32) << name << ": "
                      << std::fixed << std::setw(12) << std::setprecision(8)
                      << timers_map[name].get_duration();
                if (reference_timer != "")
                {
                     std::cout << std::right << std::setw(10) << std::setprecision(2) << 100*timers_map[name].get_duration()/timers_map[reference_timer].get_duration() << " \%"<< std::endl;
                }
                else
                {
                    std::cout << std::endl;
                }
            }
            if (reference_timer != "")
            {
                std::cout << std::string(57, '-') << std::endl;
            }
            else
            {
                std::cout << std::string(46, '-') << std::endl;
            }    
        }

        /**
         * @brief prints a comma-separated list table with durations recorded by all requested timers. This is done without any neat formatting for output post-processing from the log of the programme.
         * @param timers_names the names of the timers which durations are to be printed.
         */
        void log_timers(std::vector<std::string> timers_names)
        {
            std::cout << std::setprecision(8);
            // Read timers names requested
            for (int i; i < timers_names.size(); ++i)
            {
                if (i == (timers_names.size() - 1))
                {
                    std::cout << timers_names[i] << std::endl;
                } else {
                    std::cout << timers_names[i] << ",";    
                }
            } 

            // Read the recorded times
            for (int i; i < timers_names.size(); ++i)
            {
                if (i == (timers_names.size() - 1))
                {
                    std::cout << timers_map[timers_names[i]].get_duration() << std::endl;
                } else {
                    std::cout << timers_map[timers_names[i]].get_duration() << ",";
                }
            }
        }
        #ifdef WITH_MPI
        /**
         * @brief collects all the time data from all ranks on master rank 0, thus populating the \ref rank_durations_map object.
         * @param timers_names the names of the timers which durations are to be printed.
         */
        void collect_parallel_timers(std::vector<std::string> timers_names)
        {
            int num_timers = timers_names.size();
            std::vector<double> durations_vector(num_timers*num_ranks, 0.0);

            std::vector<double> collected_durations_vector;
            if (rank == 0)
            {
                collected_durations_vector.resize(num_timers*num_ranks, 0.0);
            }

            // Write the recorded durations into the duration_vector
            for (int timer_i = 0; timer_i < num_timers; ++timer_i)
            {
                durations_vector[rank*num_timers + timer_i] = timers_map[timers_names[timer_i]].get_duration();
            }

            double* timer_data_start = durations_vector.data()+num_timers*rank;
            MPI_Gather(timer_data_start, num_timers, MPI_DOUBLE, collected_durations_vector.data(), num_timers, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (rank == 0) 
            {
            rank_durations_map.clear();
            // fill out the rank_durations_map with content from the gatehred durations_vector
            for (int rank_i = 0; rank_i <num_ranks; ++rank_i)
            {
              rank_durations_map[rank_i] = std::vector<double>(collected_durations_vector.begin()+rank_i*num_timers, collected_durations_vector.begin()+(rank_i + 1)*num_timers);   
            }    
            }
        }

        /**
         * @brief Collects the timers on one master rank 0, then prints a comma-separated list table with durations recorded by all requested timers from all ranks.
         * @param timers_names the names of the timers which durations are to be printed.
         */
        void log_parallel_timers(std::vector<std::string>& timers_names)
        {
            collect_parallel_timers(timers_names);
            if (rank == 0)
            {
                read_rank_durations_map(timers_names);
            }
            
        }
        #endif 
        
        /**
         * @brief prints the \ref rank_durations_map object.
         * @param timers_names the names of the timers which durations are to be printed.
         */
        void read_rank_durations_map(std::vector<std::string> timers_names)
        {
            int num_timers = timers_names.size();
            // prepare the output to include the rank as well.
            timers_names.insert(timers_names.begin(), "rank");
            std::cout << std::setprecision(8);
            // Read timers names requested
            for (int i; i < timers_names.size(); ++i)
            {
                if (i == (timers_names.size() - 1))
                {
                    std::cout << timers_names[i] << std::endl;
                } else {
                    std::cout << timers_names[i] << ",";    
                }
            } 

            for (int rank_i = 0; rank_i <num_ranks; ++rank_i)
            {
                // Read the recorded times
                for (int i = 0; i < num_timers; ++i)
                {
                    if (i == 0) 
                    {
                        std::cout << rank_i << ",";
                    }
                    if (i == (num_timers - 1))
                    {
                        std::cout << rank_durations_map[rank_i][i] << std::endl;
                    } else {
                        std::cout << rank_durations_map[rank_i][i] << ",";
                    }
                }
            }      
        }

        /**
         * @brief Get the rank durations map object
         * 
         */
        std::map<int, std::vector<double>> get_rank_durations_map() 
        {
            if (rank != 0) 
            {
                std::cout << "Only rank 0 can call TimeKeeper::get_rank_durations_map(), but was called by rank " << rank << std::endl;
                exit(1);
            }
            return rank_durations_map;
        }
};
#endif