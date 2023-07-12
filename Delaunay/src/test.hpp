#ifndef __TEST_DELAUNAY_H
#define __TEST_DELAUNAY_H

#include <gtest/gtest.h>
#include "iostream"
#include "empty_class.hpp"

TEST(TestSorting, TestMergeSort)    //DA RIVEDERE NOMI!!!!
{
  vector<int> v = {44, 25, 10, 31, 25, 48, 37, 43, 18, 48, 27};
  MergeSort<int>(v, 0, v.size()-1); 
  vector<int> sortedV = {10, 18, 25, 25, 27, 31, 37, 43, 44, 48, 48};
  EXPECT_EQ(v, sortedV);
}

TEST()


