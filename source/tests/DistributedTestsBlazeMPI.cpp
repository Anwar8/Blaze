#include "DistributedMeshTests_MPI.hpp"
#include "DistributedManagersTests_MPI.hpp"
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}
