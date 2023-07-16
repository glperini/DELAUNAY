#include "test_Point.hpp"
#include "test_Segment.hpp"
#include "test_Triangle.hpp"
#include "test_Delaunay.hpp"
#include "test_TriangularMesh.hpp"
#include "test_HeapSort.hpp"

#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
};
