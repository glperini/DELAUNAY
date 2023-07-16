#ifndef __TEST_POINT_H
#define __TEST_POINT_H

#include <gtest/gtest.h>
#include "Point.hpp"

using namespace ProjectLibrary;

using namespace testing;

///TEST POINT
TEST(TestPoint, Test_getDistance)
{
  Point p1 = Point({1.0, 1.0});
  Point p2 = Point({5.0, 4.0});

  double distance = p1.getDistance(p2);
  double expectedDistance = 5.0;

  EXPECT_EQ(distance, expectedDistance);
}

TEST(TestPoint, Test_getOrder_ByX) {

  // Creo due punti da confrontare
  Point point1 = Point({2.0, 4.0});
  Point point2 = Point({3.0, 1.0});

  bool sortingByX = true;  // Ordinamento per X

  // Uso getOrder
  Point result = point1.getOrder(point2, sortingByX);

  // Verifico che il punto restituito sia uguale a point1, perchè point1.x <= point2.x
  EXPECT_EQ(result.areSamePoint(point1), true);
}

TEST(TestPoint, Test_getOrder_ByY) {

  // Creo due punti da confrontare
  Point point1 = Point({2.0, 4.0});
  Point point2 = Point({3.0, 1.0});

  bool sortingByX = false;  // Ordinamento per X

  // Uso getOrder
  Point result = point1.getOrder(point2, sortingByX);

  // Verifico che il punto restituito sia uguale a point1, perchè point1.x <= point2.x
  EXPECT_EQ(result.areSamePoint(point2), true);
}

//Altro test per getOrder nel caso in cui x uguali, controllo se va a ordinare sulle y
TEST(TestPoint, Test_getOrder_ByXSameX) {

    // Definisco due punti con coordinate X uguali
    Point point1 = Point({2.0, 4.0});
    Point point2 = Point({2.0, 3.0});

    bool sortingByX = true; // Ordino per coordinata X

    // Uso getOrder
    Point result = point1.getOrder(point2, sortingByX);

    // Verifica che il punto restituito sia quello corretto (point2.y < point1.y)
    EXPECT_EQ(result.areSamePoint(point2), true);
}

//Altro test per getOrder nel caso in cui y uguali, controllo se va a ordinare sulle x
TEST(TestPoint, Test_getOrder_ByYSameY) {

    // Definisco due punti con coordinate X uguali
    Point point1 = Point({2.0, 4.0});
    Point point2 = Point({1.0, 4.0});

    bool sortingByX = false; // Ordino per coordinata Y

    // Uso getOrder
    Point result = point1.getOrder(point2, sortingByX);

    // Verifica che il punto restituito sia quello corretto (point2.x < point1.x)
    EXPECT_EQ(result.areSamePoint(point2), true);
}

#endif // __TEST_POINT_H
