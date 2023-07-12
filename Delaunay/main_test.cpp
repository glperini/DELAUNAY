#include "test_Sorting.hpp"
#include "test_Objects.hpp"
#include "test_Delaunay.hpp"

#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
