#ifndef TIME_KEEPERS_TESTS_HPP
#define TIME_KEEPERS_TESTS_HPP

#include "TestHelpers.hpp"

class TimeKeeperParallelTests : public ::testing::Test {
    public:
      TimeKeeper time_keeper;
      int rank, num_ranks;
      std::vector<std::string> timer_names = {"test_1", "test_2", "test_3", "test_4", "test_5"};
      
      void SetUp() override {
            get_my_rank(rank);
            get_num_ranks(num_ranks);
            time_keeper.initialise_parallel_keeper(rank, num_ranks);
            time_keeper.add_timers(timer_names);
        }
      void TearDown() override {
  }
  };
  TEST_F(TimeKeeperParallelTests, timer_gather)
  {  
    time_keeper.start_timer(timer_names[rank]);
    sleep(1);
    time_keeper.stop_timer(timer_names[rank]);
    time_keeper.log_parallel_timers(timer_names);
    
    if (rank == 0)
    {
        std::map<int, std::vector<double>> durations_map = time_keeper.get_rank_durations_map();
        for (int rank_i = 0; rank_i < num_ranks; ++rank_i)
        {
            for (int timer_i = 0; timer_i < timer_names.size(); ++timer_i)
            {
                if (timer_i == rank_i)
                {
                    EXPECT_GT(durations_map[rank_i][rank_i], 1.0);
                } else {
                    EXPECT_NEAR(durations_map[rank_i][timer_i], 0.0, 1e-8);
                }
            }
        }
    }
    
  }
#endif