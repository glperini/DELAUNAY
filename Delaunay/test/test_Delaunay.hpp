#ifndef __TEST_DELAUNAY_H
#define __TEST_DELAUNAY_H

#include <gtest/gtest.h>
#include "objects.hpp"
#include "sorting.hpp"
#include "delaunay.hpp"
#include "functions.hpp"

using namespace LibraryObjects;
using namespace LibrarySorting;
using namespace LibraryDelaunay;
using namespace LibraryFunctions;

using namespace testing;

///TEST DELAUNAY
TEST(TestDelaunay, Test_isPointInsideCircumcircle_In) {
    Delaunay delaunay;

    vector<Triangle> triangles;
    // Creo  triangoli di prova
    Point p1 = Point(0.0, 0.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(0.0, 3.0);
    Point p4 = Point(3.0, 3.0);
    Point p5 = Point(5.0, 0.0);

    triangles = {Triangle(p1, p2, p3), Triangle(p2, p4, p3), Triangle(p2, p5, p4)};

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideCircle = Point(4.0, 1.0);

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideCircumcircle(pointInsideCircle, triangles, triangleIndex);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, 2);
}

TEST(TestDelaunay, Test_isPointInsideCircumcircle_Out) {
    Delaunay delaunay;

    vector<Triangle> triangles;
    // Creo  triangoli di prova
    Point p1 = Point(0.0, 0.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(0.0, 3.0);
    Point p4 = Point(3.0, 3.0);
    Point p5 = Point(5.0, 0.0);

    triangles = {Triangle(p1, p2, p3), Triangle(p2, p4, p3), Triangle(p2, p5, p4)};

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideCircle = Point(10.0, 0.0);

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideCircumcircle(pointInsideCircle, triangles, triangleIndex);

    // Verifico che la funzione abbia restituito un valore diverso da -1
    EXPECT_EQ(resultIndex, -1);
}


TEST(TestDelaunay, Test_isPointInsideTriangle_In) {
    Delaunay delaunay;

    vector<Triangle> triangles;
    // Creo  triangoli di prova
    Point p1 = Point(0.0, 0.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(0.0, 3.0);

    triangles.push_back(Triangle(p1, p2, p3));

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideTriangle = Point(1.0, 1.0);

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideTriangle(pointInsideTriangle, triangles, triangleIndex);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, 0);
}

TEST(TestDelaunay, Test_isPointInsideTriangle_Out) {

    Delaunay delaunay;

    vector<Triangle> triangles;
    // Creo  triangoli di prova
    Point p1 = Point(0.0, 0.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(0.0, 3.0);

    triangles.push_back(Triangle(p1, p2, p3));

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideTriangle = Point(2.0, 2.0);

    //chiamo la funzione
   int resultIndex = delaunay.isPointInsideTriangle(pointInsideTriangle, triangles, triangleIndex);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, -1);
}

TEST(TestDelaunay, Test_findTriangleContainingPoint_In) {

    Delaunay delaunay;

    vector<Triangle> triangles;
    //Creo triangoli di prova
    Point p1 = Point(2.0, 1.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(5.0, 1.5);
    Point p4 = Point(4.0, 2.0);
    Point p5 = Point(5.0, 0.0);

    triangles = {Triangle(p1, p2, p3), Triangle(p2, p5, p3), Triangle(p1, p3, p4)};

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointOutsideTriangle = Point(4.0, 1.5);

    //chiamo la funzione
    int resultIndex = delaunay.findTriangleContainingPoint(pointOutsideTriangle, triangles, triangleIndex);

    EXPECT_EQ(resultIndex, 2);
}

TEST(TestDelaunay, Test_findTriangleContainingPoint_Out) {

    Delaunay delaunay;

    vector<Triangle> triangles;
    //Creo triangoli di prova
    Point p1 = Point(2.0, 1.0);
    Point p2 = Point(3.0, 0.0);
    Point p3 = Point(5.0, 1.5);
    Point p4 = Point(4.0, 2.0);
    Point p5 = Point(5.0, 0.0);

    triangles = {Triangle(p1, p2, p3), Triangle(p2, p5, p3), Triangle(p1, p3, p4)};

    int triangleIndex = 0;
    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointOutsideTriangle = Point(2.0, 0.0);

    //chiamo la funzione
    int resultIndex = delaunay.findTriangleContainingPoint(pointOutsideTriangle, triangles, triangleIndex);

    EXPECT_EQ(resultIndex, -1);
}

//verifichiamo caso: adjacentSide1 = 1, adjacentSide2 = 2
TEST(TestDelaunay, Test_flipTriangles_Side12Neighbors) {

    Delaunay delaunay;
    Point p1 = Point(1.0, 2.0);
    Point p2 = Point(2.0, 0.0);
    Point p3 = Point(2.0, 4.0);
    Point p4 = Point(3.5, 2.0);
    Point p5 = Point(1.0, 0.0);
    Point p6 = Point(3.5, 4.5);

    Triangle triangle1 = Triangle(p1,p2,p3);
    Triangle triangle2 = Triangle(p2,p4,p3);
    Triangle triangle3 = Triangle(p1,p5,p2);
    Triangle triangle4 = Triangle(p4,p6,p3);
    vector<Triangle> triangles = {triangle1, triangle2, triangle3, triangle4};

    int adjacentSide1 = 1;
    int adjacentSide2 = 2;
    vector<int> hullTrianglesIndices = {0, 1, 2, 3};
    int triangle1Index = 0;
    int triangle2Index = 1;
    triangles[0].neighbors = {2,1,-1};
    triangles[1].neighbors = {-1,3,0};
    triangles[2].neighbors = {-1,-1,0};
    triangles[3].neighbors = {-1,-1,1};


    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2, hullTrianglesIndices, triangles);

    vector<int> expectedNeighbors1 = {2,-1,1};
    vector<int> expectedNeighbors2 = {3,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedNeighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedNeighbors2[i]);
    }
}

//verifichiamo altro caso: adjacentSide1 = 2, adjacentSide2 = 1
TEST(TestDelaunay, Test_flipTriangles_Side21Neighbors) {

    Delaunay delaunay;
    Point p1 = Point(1.0, 2.0);
    Point p2 = Point(2.0, 0.0);
    Point p3 = Point(2.0, 4.0);
    Point p4 = Point(3.5, 2.0);
    Point p5 = Point(1.0, 0.0);
    Point p6 = Point(3.5, 4.5);

    Triangle triangle1 = Triangle(p1,p2,p4);
    Triangle triangle2 = Triangle(p3,p1,p4);
    Triangle triangle3 = Triangle(p1,p5,p2);
    Triangle triangle4 = Triangle(p4,p6,p3);
    vector<Triangle> triangles = {triangle1, triangle2, triangle3, triangle4};

    int adjacentSide1 = 2;
    int adjacentSide2 = 1;
    vector<int> hullTrianglesIndices = {0, 1, 2, 3};
    int triangle1Index = 0;
    int triangle2Index = 1;
    triangles[0].neighbors = {2,-1,1};
    triangles[1].neighbors = {-1,0,3};
    triangles[2].neighbors = {-1,-1,0};
    triangles[3].neighbors = {-1,-1,1};

    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2, hullTrianglesIndices, triangles);

    vector<int> expectedNeighbors1 = {-1,3,1};
    vector<int> expectedNeighbors2 = {-1,2,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedNeighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedNeighbors2[i]);
    }
}

//verifichiamo l'aggiornamento di hullTriangleIndices
TEST(TestDelaunay, Test_flipTriangles_hullTriangleIndices) {

    Delaunay delaunay;
    Point p1 = Point(1.0, 2.0);
    Point p2 = Point(2.0, 0.0);
    Point p3 = Point(2.0, 4.0);
    Point p4 = Point(3.5, 2.0);
    Point p5 = Point(1.0, 0.0);
    Point p6 = Point(3.5, 4.5);
    Point p7 = Point(3.0, 0.0);

    Triangle triangle1 = Triangle(p1,p2,p3);
    Triangle triangle2 = Triangle(p2,p4,p3);
    Triangle triangle3 = Triangle(p1,p5,p2);
    Triangle triangle4 = Triangle(p4,p6,p3);
    Triangle triangle5 = Triangle(p2,p7,p4);
    vector<Triangle> triangles = {triangle1, triangle2, triangle3, triangle4, triangle5};

    int adjacentSide1 = 1;
    int adjacentSide2 = 2;
    vector<int> hullTrianglesIndices = {0, 2, 3, 4};
    int triangle1Index = 0;
    int triangle2Index = 1;
    triangles[0].neighbors = {2,1,-1};
    triangles[1].neighbors = {4,3,0};
    triangles[2].neighbors = {-1,-1,0};
    triangles[3].neighbors = {-1,-1,1};
    triangles[4].neighbors = {-1,-1,1};

    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2, hullTrianglesIndices, triangles);

    vector<int> expectedHullIndices = {2,3,4,1};
    for (int i=0;i<4;i++)
        EXPECT_EQ(hullTrianglesIndices[i],expectedHullIndices[i]);
}

//
TEST(TestDelaunay, Test_verifyDelaunayCondition) {

    Delaunay delaunay;
    Point p1 = Point(1.0, 2.0);
    Point p2 = Point(2.0, 0.0);
    Point p3 = Point(2.0, 4.0);
    Point p4 = Point(3.5, 2.0);
    Point p5 = Point(3.5, 4.5);
    Point p6 = Point(3.0, 0.0);

    Triangle triangle1 = Triangle(p1,p2,p3);
    Triangle triangle2 = Triangle(p2,p4,p3);
    Triangle triangle3 = Triangle(p4,p5,p3);
    Triangle triangle4 = Triangle(p2,p6,p4);
    vector<Triangle> triangles = {triangle1, triangle2, triangle3, triangle4};

    vector<int> hullTriangleIndices = {0,2,3};

    int triangleIndex = 1;
    triangles[0].neighbors = {-1,1,-1};
    triangles[1].neighbors = {3,2,0};
    triangles[2].neighbors = {-1,-1,1};
    triangles[3].neighbors = {-1,-1,1};

    delaunay.verifyDelaunayCondition(triangles[1], triangleIndex, triangles, hullTriangleIndices);

    vector<int> expectedNeighbors1 = {-1,3,1};
    vector<int> expectedNeighbors2 = {2,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedNeighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedNeighbors2[i]);
    }

    vector<int> expectedHullTriangleIndices = {0,2,3,1};
    for (int i=0; i<4; i++)
        EXPECT_EQ(expectedHullTriangleIndices[i], hullTriangleIndices[i]);
}

TEST(TestDelaunay, Test_inTriangle) {
    Delaunay delaunay;
    Point p1 = Point(3.0, 3.0);
    Point p2 = Point(5.0, 7.0);
    Point p3 = Point(11.0, 3.0);
    Point p4 = Point(6.0, 2.0);
    Point p5 = Point(3.0, 6.0);
    Point InsidePoint = Point(6.0, 5.0);

    Triangle triangle1 = Triangle(p1,p4,p2);
    Triangle triangle2 = Triangle(p2,p4,p3);
    Triangle triangle3 = Triangle(p1,p2,p5);
    vector<Triangle> triangles = {triangle1, triangle2, triangle3};

    vector<int> hullTriangleIndices = {0, 1, 2};

    int triangleIndex = 1;
    triangles[0].neighbors = {-1,1,2};
    triangles[1].neighbors = {0,-1,-1};
    triangles[2].neighbors = {0,-1,-1};

    delaunay.inTriangle(InsidePoint, triangles, triangleIndex, hullTriangleIndices);

    vector<int> expectedNeighbors1 = {-1,3,2};
    vector<int> expectedNeighbors2 = {4,-1,2};
    vector<int> expectedNeighbors3 = {-1,0,1};
    vector<int> expectedNeighbors4 = {-1,4,0};
    vector<int> expectedNeighbors5 = {-1,1,3};


    for (int i=0; i<3; i++) {
        ASSERT_EQ(triangles[0].neighbors[i], expectedNeighbors1[i]);
        ASSERT_EQ(triangles[1].neighbors[i], expectedNeighbors2[i]);
        ASSERT_EQ(triangles[2].neighbors[i], expectedNeighbors3[i]);
        ASSERT_EQ(triangles[3].neighbors[i], expectedNeighbors4[i]);
        ASSERT_EQ(triangles[4].neighbors[i], expectedNeighbors5[i]);
    }

}

TEST(TestDelaunay, Test_outTriangle) {

      Delaunay delaunay;

      vector<Triangle> triangles;
      vector<Segment> convexHull;
      vector<Point> hullPoints;
      vector<int> hullTrianglesIndices;

      Point outsidePoint = Point(-2.0, 2.0);

      Point p1 = Point(0.0, 0.0);
      Point p2 = Point(3.0, 0.0);
      Point p3 = Point(1.0, 2.0);
      Point p4 = Point(4.0, 2.0);

      triangles = {Triangle(p1, p2, p3), Triangle(p2, p4, p3)};
      triangles[0].neighbors = {-1,1,-1};
      triangles[1].neighbors = {-1,-1,0};

      Vector2d linecoeff1 = {0.0, 0.0};
      Vector2d linecoeff2 = {2.0, -6.0};
      Vector2d linecoeff3 = {0.0, 2.0};
      Vector2d linecoeff4 = {2.0, 0.0};
      convexHull = {Segment(p1, p2, linecoeff1), Segment(p2, p4, linecoeff2), Segment(p4, p3, linecoeff3), Segment(p3, p1, linecoeff4)};

      hullPoints = {p1, p2, p3, p4};

      hullTrianglesIndices = {0, 1};

      delaunay.outTriangle(outsidePoint, triangles, convexHull, hullPoints, hullTrianglesIndices);

      //verifico che il triangolo creato sia corretto
      EXPECT_EQ(areSamePoint(triangles[2].p1, outsidePoint), true);

      EXPECT_EQ(areSamePoint(triangles[2].p2, p1), true);

      EXPECT_EQ(areSamePoint(triangles[2].p3, p3), true);

      vector<int> expectedNeighbors1 = {-1, 1, 2};
      vector<int> expectedNeighbors2 = {-1, -1, 0};
      vector<int> expectedNeighbors3 = {-1, 0, -1};

      for (int i=0; i<3; i++) {
          EXPECT_EQ(triangles[2].neighbors[i], expectedNeighbors3[i]);
          EXPECT_EQ(triangles[0].neighbors[i], expectedNeighbors1[i]);
          EXPECT_EQ(triangles[1].neighbors[i], expectedNeighbors2[i]);
          }

      //verifico i punti nel guscio: expectedHullPoints = {p1, p2, p3, p4, outsidePoint};
      EXPECT_EQ(areSamePoint(hullPoints[4],outsidePoint), true);

      //verifico che abbia aggiunto il nuovo triangolo nel guscio: expectedHullTrianglesIndices = {0, 1, 2}
      EXPECT_EQ(hullTrianglesIndices[2], 2);
}


TEST(TestDelaunay, Test_drawSegments) {
    Delaunay delaunay;

    Point p1 = Point(3.0, 3.0);
    Point p2 = Point(5.0, 7.0);
    Point p3 = Point(11.0, 3.0);
    Point p4 = Point(6.0, 2.0);
    Point p5 = Point(3.0, 6.0);

    Triangle triangle1 = Triangle(p1,p3,p2);
    Triangle triangle2 = Triangle(p1,p4,p3);
    Triangle triangle3 = Triangle(p2,p5,p1);

    vector<Triangle> triangles = {triangle1, triangle2, triangle3};

    vector<Segment> segments = delaunay.drawSegments(triangles);

    vector<Segment> expectedSegments = {Segment(p1,p3), Segment(p3,p2), Segment(p2,p1), Segment(p1,p4), Segment(p4,p3), Segment(p2,p5), Segment(p5,p1)};

    for (int i = 0; i<7; i++)
        EXPECT_EQ(areSameSegment(segments[i], expectedSegments[i]), true);
}

#endif // __TEST_DELAUNAY_H
