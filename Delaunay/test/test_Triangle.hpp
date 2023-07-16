#ifndef __TEST_TRIANGLE_H
#define __TEST_TRIANGLE_H

#include <gtest/gtest.h>
#include "Triangle.hpp"

using namespace ProjectLibrary;

using namespace testing;

///TEST TRIANGLE
TEST(TestTriangle, Test_antiClockWiseOrder)  {

    //Definisco tre punti non ordinati in senso antiorario
    Point p1 = Point({0.0, 0.0});
    Point p2 = Point({1.0, 1.0});
    Point p3 = Point({2.0, 0.0});

    // Creo Triangle
    Triangle triangle = Triangle({p1, p2, p3});

    // Chiamo funzione anticlockwise
    triangle.antiClockWiseOrder();

    // Verifico che i punti siano in senso antiorario
    Point expectedP1 = Point({0.0, 0.0});  //punti giusti in senso antiorario
    Point expectedP2 = Point({2.0, 0.0});
    Point expectedP3 = Point({1.0, 1.0});

    EXPECT_EQ(triangle.get_point(0).areSamePoint(expectedP1), true);
    EXPECT_EQ(triangle.get_point(1).areSamePoint(expectedP2), true);
    EXPECT_EQ(triangle.get_point(2).areSamePoint(expectedP3), true);
}


#endif // __TEST_TRIANGLE_H
