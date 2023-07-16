#ifndef __TEST_HEAPSORT_H
#define __TEST_HEAPSORT_H

#include <gtest/gtest.h>
#include "Delaunay.hpp"
#include "TriangularMesh.hpp"
#include "HeapSort.hpp"

using namespace ProjectLibrary;
using namespace SortingLibrary;

using namespace testing;

///TEST HeapSort
TEST(TestSorting, Test_heapSort_ByXAllDistinct) {

  // vettore di punti non ordinato
  vector<Point> points = {Point({2.0, 4.0}), Point({1.0, 3.0}), Point({5.0, 2.0}), Point({3.0, 1.0})};

  bool sortingByX = true;  // Ordino per X

  // funzione HeapSort
  points = heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point({1.0, 3.0}), Point({2.0, 4.0}), Point({3.0, 1.0}), Point({5.0, 2.0})};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(points[i].areSamePoint(expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByYAllDistinct) {

  // vettore di punti non ordinato
  vector<Point> points = {Point({2.0, 4.0}), Point({1.0, 3.0}), Point({5.0, 2.0}), Point({3.0, 1.0})};

  bool sortingByX = false;  // Ordino per Y

  // funzione HeapSort
  points = heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point({3.0, 1.0}), Point({5.0, 2.0}), Point({1.0, 3.0}), Point({2.0, 4.0})};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(points[i].areSamePoint(expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByXSameX) {

  // vettore di punti non ordinato
  vector<Point> points = {Point({2.0, 4.0}), Point({1.0, 3.0}), Point({2.0, 2.0}), Point({3.0, 1.0})};

  bool sortingByX = true;  // Ordino per X

  // funzione HeapSort
  points = heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point({1.0, 3.0}), Point({2.0, 2.0}), Point({2.0, 4.0}), Point({3.0, 1.0})};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(points[i].areSamePoint(expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByYSameY) {

  // vettore di punti non ordinato
  vector<Point> points = {Point({2.0, 4.0}), Point({1.0, 3.0}), Point({5.0, 3.0}), Point({3.0, 1.0})};

  bool sortingByX = false;  // Ordino per Y

  // funzione HeapSort
  points = heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point({3.0, 1.0}), Point({1.0, 3.0}), Point({5.0, 3.0}), Point({2.0, 4.0})};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(points[i].areSamePoint(expected[i]), true);
  }
}


#endif // __TEST_HEAPSORT_H
