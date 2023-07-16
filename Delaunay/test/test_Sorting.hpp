#ifndef __TEST_SORTING_H
#define __TEST_SORTING_H

#include <gtest/gtest.h>
#include "objects.hpp"
#include "sorting.hpp"
#include "functions.hpp"

using namespace LibraryObjects;
using namespace LibrarySorting;
using namespace LibraryFunctions;

using namespace testing;


///TEST SORTING
TEST(TestSorting, Test_heapSort_ByXAllDistinct) {
  HeapSorter sorter;

  // vettore di punti non ordinato
  vector<Point> points = {Point(2.0, 4.0), Point(1.0, 3.0), Point(5.0, 2.0), Point(3.0, 1.0)};

  bool sortingByX = true;  // Ordino per X

  // funzione HeapSort
  points = sorter.heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point(1.0, 3.0), Point(2.0, 4.0), Point(3.0, 1.0), Point(5.0, 2.0)};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(areSamePoint(points[i], expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByYAllDistinct) {
  HeapSorter sorter;

  // vettore di punti non ordinato
  vector<Point> points = {Point(2.0, 4.0), Point(1.0, 3.0), Point(5.0, 2.0), Point(3.0, 1.0)};

  bool sortingByX = false;  // Ordino per Y

  // funzione HeapSort
  points = sorter.heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point(3.0, 1.0), Point(5.0, 2.0), Point(1.0, 3.0), Point(2.0, 4.0)};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(areSamePoint(points[i], expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByXSameX) {
  HeapSorter sorter;

  // vettore di punti non ordinato
  vector<Point> points = {Point(2.0, 4.0), Point(1.0, 3.0), Point(2.0, 2.0), Point(3.0, 1.0)};

  bool sortingByX = true;  // Ordino per X

  // funzione HeapSort
  points = sorter.heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point(1.0, 3.0), Point(2.0, 2.0), Point(2.0, 4.0), Point(3.0, 1.0)};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(areSamePoint(points[i], expected[i]), true);
  }
}

TEST(TestSorting, Test_heapSort_ByYSameY) {
  HeapSorter sorter;

  // vettore di punti non ordinato
  vector<Point> points = {Point(2.0, 4.0), Point(1.0, 3.0), Point(5.0, 3.0), Point(3.0, 1.0)};

  bool sortingByX = false;  // Ordino per Y

  // funzione HeapSort
  points = sorter.heapSort(points, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point(3.0, 1.0), Point(1.0, 3.0), Point(5.0, 3.0), Point(2.0, 4.0)};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(areSamePoint(points[i], expected[i]), true);
  }
}

TEST(TestSorting, Test_getOrder_ByX) {
  HeapSorter sorter;  // Creo HeapSorter

  // Creo due punti da confrontare
  Point point1 = Point(2.0, 4.0);
  Point point2 = Point(3.0, 1.0);

  bool sortingByX = true;  // Ordinamento per X

  // Uso getOrder
  Point result = sorter.getOrder(point1, point2, sortingByX);

  // Verifico che il punto restituito sia uguale a point1, perchè point1.x <= point2.x
  EXPECT_EQ(areSamePoint(result, point1), true);
}

TEST(TestSorting, Test_getOrder_ByY) {
  HeapSorter sorter;  // Creo HeapSorter

  // Creo due punti da confrontare
  Point point1 = Point(2.0, 4.0);
  Point point2 = Point(3.0, 1.0);

  bool sortingByX = false;  // Ordinamento per X

  // Uso getOrder
  Point result = sorter.getOrder(point1, point2, sortingByX);

  // Verifico che il punto restituito sia uguale a point1, perchè point1.x <= point2.x
  EXPECT_EQ(areSamePoint(result, point2), true);
}

//Altro test per getOrder nel caso in cui x uguali, controllo se va a ordinare sulle y
TEST(TestSorting, Test_getOrder_ByXSameX) {
    HeapSorter sorter; // Creo HeapSorter

    // Definisco due punti con coordinate X uguali
    Point point1 = Point(2.0, 4.0);
    Point point2 = Point(2.0, 3.0);

    bool sortingByX = true; // Ordino per coordinata X

    // Uso getOrder
    Point result = sorter.getOrder(point1, point2, sortingByX);

    // Verifica che il punto restituito sia quello corretto (point2.y < point1.y)
    EXPECT_EQ(areSamePoint(result, point2), true);
}

//Altro test per getOrder nel caso in cui y uguali, controllo se va a ordinare sulle x
TEST(TestSorting, Test_getOrder_ByYSameY) {
    HeapSorter sorter; // Creo HeapSorter

    // Definisco due punti con coordinate X uguali
    Point point1 = Point(2.0, 4.0);
    Point point2 = Point(1.0, 4.0);

    bool sortingByX = false; // Ordino per coordinata Y

    // Uso getOrder
    Point result = sorter.getOrder(point1, point2, sortingByX);

    // Verifica che il punto restituito sia quello corretto (point2.x < point1.x)
    EXPECT_EQ(areSamePoint(result, point2), true);
}

#endif // __TEST_SORTING_H
