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
    public:
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
            
            std::cout << std::endl <<  "Reading some ExecutionTimers:" << std::endl;
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

};
#endif