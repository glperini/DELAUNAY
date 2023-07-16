#ifndef __TEST_SEGMENT_H
#define __TEST_SEGMENT_H

#include <gtest/gtest.h>
#include "Segment.hpp"

using namespace ProjectLibrary;

using namespace testing;

///TEST SEGMENT
TEST(TestSegment, Test_lineCoefficients_SameX) {
    Point p1 = Point({2.0, 3.0});
    Point p2 = Point({2.0, 4.0});
    Segment segment = Segment({p1, p2});
    segment.lineCoefficients();

    EXPECT_TRUE(isnan(segment.get_coefficient(0)));       // Verifica se coefficients[0] è NaN
    EXPECT_EQ(segment.get_coefficient(1), 2.0);           // Verifica se coefficients[1] è uguale a p1.x
}

TEST(TestSegment , Test_lineCoefficients) {
    Point p1 = Point({1.0, 2.0});
    Point p2 = Point({3.0, 4.0});
    Segment segment = Segment({p1, p2});
    segment.lineCoefficients();
    EXPECT_DOUBLE_EQ(segment.get_coefficient(0), 1.0);    // Verifica se coefficients[0] è uguale a 1.0
    EXPECT_DOUBLE_EQ(segment.get_coefficient(1), 1.0);    // Verifica se coefficients[1] è uguale a 1.0
}

TEST(TestSegment, Test_getAngle_ParallelSegments)
{
    //Costruiamo due coppie di punti, per poter costruire i due segmenti
    Point s1 = Point({1.0, 1.0});  //Primo punto, primo segmento
    Point s2 = Point({1.0, 5.0});  //Secondo punto, primo segmento
    Point t1 = Point({3.0, 0.0});  //Primo punto, secondo segmento
    Point t2 = Point({3.0, 2.0});  //Secondo punto, secondo segmento

    //Costruiamo i segmenti
    Segment segment_S = Segment({s1, s2});
    Segment segment_T = Segment({t1, t2});

    //Calcoliamo i coefficienti delle rette associate ai due segmenti
    segment_S.lineCoefficients();
    segment_T.lineCoefficients();

    //Confrontiamo gli angoli
    double expectedAngleDegree = 0.0;
    double angleDegree = segment_S.getAngle(segment_T);

    EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}

TEST(TestSegment, Test_getAngle_PerpendicularSegments)
{
    //Costruiamo due coppie di punti, per poter costruire i due segmenti
    Point s1 = Point({1.0, 1.0});  //Primo punto, primo segmento
    Point s2 = Point({3.0, 3.0});  //Secondo punto, primo segmento
    Point t1 = Point({4.0, 5.0});  //Primo punto, secondo segmento
    Point t2 = Point({6.0, 3.0});  //Secondo punto, secondo segmento

    //Costruiamo i segmenti
    Segment segment_S = Segment({s1, s2});
    Segment segment_T = Segment({t1, t2});

    //Calcoliamo i coefficienti delle rette associate ai due segmenti
    segment_S.lineCoefficients();
    segment_T.lineCoefficients();

    //Confrontiamo gli angoli
    double expectedAngleDegree = 90.0;
    double angleDegree = segment_S.getAngle(segment_T);

    EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}

TEST(TestSegment, Test_getAngle_GenericSegments)
{
    //Costruiamo due coppie di punti, per poter costruire i due segmenti
    Point s1 = Point({1.0, 1.0});  //Primo punto, primo segmento
    Point s2 = Point({4.0, 1.0});  //Secondo punto, primo segmento
    Point t1 = Point({1.0, 2.0});  //Primo punto, secondo segmento
    Point t2 = Point({3.0, 4.0});  //Secondo punto, secondo segmento


    //Costruiamo i segmenti
    Segment segment_S = Segment({s1, s2});
    Segment segment_T = Segment({t1, t2});

    //Calcoliamo i coefficienti delle rette associate ai due segmenti
    segment_S.lineCoefficients();
    segment_T.lineCoefficients();

    //Confrontiamo gli angoli
    double expectedAngleDegree = 45.0;
    double angleDegree = segment_S.getAngle(segment_T);

    EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}

#endif // __TEST_SEGMENT_H
