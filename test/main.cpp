/**
 * @authors     Remo Diethelm
 * @affiliation ANYbotics
 * @brief       Main test file.
 */

#include <gtest/gtest.h>

// Run all the tests that were declared with TEST()
int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  const unsigned int currentTime = time(nullptr);
  std::cerr << "Seed for random numbers: " << currentTime << std::endl;  // cerr output is printed in GTests.
  srand(currentTime);
  return RUN_ALL_TESTS();
}
