/**
 * @file ExecutionTimer.hpp
 * @brief Defines the ExecutionTimer class which measures the duration of code execution and used for profiling.
 */
#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER
#include <sys/time.h>
#include <string>
#include "maths_defaults.hpp"

/**
 * @brief Times execution of the code. 
 * @details Correct usage of the timer is to call \ref start and \ref stop, and the timer will automatically measure the duration between the two. The \ref stop command does not terminate the timer, and the timer can be started and stopped again with the duration between each start-stop pair added to the timer duration. To reset timer to zero, use the \ref reset function.
 */
class ExecutionTimer 
{
    protected:
        double duration = 0.0; /**< cumulative duration between all start-stop pairs.*/
        double start_time = -1.0; /**< time at which the \ref start function was last called. */
        double stop_time = -1.0; /**< time at which the \ref stop function was last called. */
    public:
        /**
         * @brief Get the current time using the `sys/time.h` header.
         * 
         * @return double current time.
         */
        double get_current_time()
        {
            struct timeval tp;
            gettimeofday (&tp, NULL);
            return tp.tv_sec + tp.tv_usec/(double)1.0e6;
        }

        /**
         * @brief set the \ref start_time to current time using \ref get_current_time.
         * 
         */
        void start()
        {
            start_time = get_current_time();
        }

        /**
         * @brief set the \ref stop_time to current time using \ref get_current_time.
         * 
         */
        void stop()
        {
            stop_time = get_current_time();
            duration += stop_time - start_time;
        }

        /**
         * @brief resets the timer setting \ref duration to 0.0, and \ref start_time and \ref stop_time to -1.0.
         * 
         */
        void reset()
        {
            start_time = stop_time = -1.0;
            duration = 0.0; 
        }

        /**
         * @brief Get the duration value this timer has mesured so far.
         * 
         * @return double \ref duration.
         */
        double get_duration() const {return duration;}
};
#endif 