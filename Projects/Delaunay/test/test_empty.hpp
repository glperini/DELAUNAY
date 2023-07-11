#ifndef __TEST_EMPTY_H
#define __TEST_EMPTY_H

#include <gtest/gtest.h>
#include "empty_class.hpp"

using namespace testing;
using namespace ProjectLibrary;

///TEST SORTING
TEST(TestSorting, Test_MergeSort) {
  MergeSorter sorter;

  // vettore di punti non ordinato
  vector<Point> points = {Point(2.0, 4.0), Point(1.0, 3.0), Point(5.0, 2.0), Point(3.0, 1.0)};

  bool sortingByX = true;  // Ordino per X

  // funzione MergeSort
  int sx = 0;
  int dx = points.size() - 1;
  points = sorter.MergeSort(points, sx, dx, sortingByX);

  // Verifica che il vettore sia stato ordinato giusto
  vector<Point> expected = {Point(1.0, 3.0), Point(2.0, 4.0), Point(3.0, 1.0), Point(5.0, 2.0)};
  for (unsigned int i=0; i<points.size(); i++) {
      EXPECT_EQ(areSamePoint(points[i], expected[i]), true);
  }
}

TEST(TestSorting, Test_getOrder) {
  MergeSorter sorter;  // Creo MergeSorter

  // Creo due punti da confrontare
  Point point1 = Point(2.0, 4.0);
  Point point2 = Point(3.0, 1.0);

  bool sortingByX = true;  // Ordinamento per X

  // Uso getOrder
  Point result = sorter.getOrder(point1, point2, sortingByX);

  // Verifico che il punto restituito sia uguale a point1, perchè point1.x <= point2.x
  EXPECT_EQ(areSamePoint(result, point1), true);
}


//Altro test per getOrder nel caso in cui x uguali, controllo se va a ordinare sulle y
TEST(TestSorting, Test_getOrder_SameX) {
    MergeSorter sorter; // Creo MergeSorter

    // Definisco due punti con coordinate X uguali
    Point point1 = Point(2.0, 4.0);
    Point point2 = Point(2.0, 3.0);

    bool sortingByX = true; // Ordino per coordinata X

    // Uso getOrder
    Point result = sorter.getOrder(point1, point2, sortingByX);

    // Verifica che il punto restituito sia quello corretto (point2.y < point1.y)
    EXPECT_EQ(areSamePoint(result, point2), true);
}

///TEST POINT
TEST(TestPoint, TestGetDistance)
{
  Point p1 = Point(1.0, 1.0);
  Point p2 = Point(5.0, 4.0);

  double distance = getDistance(p1, p2);
  double expectedDistance = 5.0;

  EXPECT_EQ(distance, expectedDistance);
}


///TEST SEGMENT
TEST(TestSegment, Test_LineCoefficients_SameX) {
    Segment segment;
    Point p1 = Point(2.0, 3.0);
    Point p2 = Point(2.0, 4.0);
    Vector2d result = segment.lineCoefficients(p1, p2);
    EXPECT_TRUE(isnan(result[0]));       // Verifica se coefficients[0] è NaN
    EXPECT_EQ(result[1], 2.0);           // Verifica se coefficients[1] è uguale a p1.x
}

TEST(TestSegment , Test_LineCoefficients) {
    Segment segment;
    Point p1 = Point(1.0, 2.0);
    Point p2 = Point(3.0, 4.0);
    Vector2d result = segment.lineCoefficients(p1, p2);
    EXPECT_DOUBLE_EQ(result[0], 1.0);    // Verifica se coefficients[0] è uguale a 1.0
    EXPECT_DOUBLE_EQ(result[1], 1.0);    // Verifica se coefficients[1] è uguale a 1.0
}

TEST(TestSegment, TestGetAngleWithParallelSegments)
{
  Segment segment;

  //Costruiamo due coppie di punti, per poter costruire i due segmenti
  Point s1 = Point(1.0, 1.0);  //Primo punto, primo segmento
  Point s2 = Point(1.0, 5.0);  //Secondo punto, primo segmento
  Point t1 = Point(3.0, 0.0);  //Primo punto, secondo segmento
  Point t2 = Point(3.0, 2.0);  //Secondo punto, secondo segmento

  //Calcoliamo i coefficienti delle rette associate ai due segmenti
  Vector2d coefficients_S = segment.lineCoefficients(s1, s2);
  Vector2d coefficients_T = segment.lineCoefficients(t1, t2);

  //Costruiamo i segmenti
  Segment segment_S = Segment(s1, s2, coefficients_S);
  Segment segment_T = Segment(t1, t2, coefficients_T);

  //Confrontiamo gli angoli
  double expectedAngleDegree = 0.0;
  double angleDegree = getAngle(segment_S, segment_T);
  EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}

TEST(TestSegment, TestGetAngleWithPerpendicularSegments)
{
  Segment segment;

  //Costruiamo due coppie di punti, per poter costruire i due segmenti
  Point s1 = Point(1.0, 1.0);  //Primo punto, primo segmento
  Point s2 = Point(3.0, 3.0);  //Secondo punto, primo segmento
  Point t1 = Point(4.0, 5.0);  //Primo punto, secondo segmento
  Point t2 = Point(6.0, 3.0);  //Secondo punto, secondo segmento

  //Calcoliamo i coefficienti delle rette associate ai due segmenti
  Vector2d coefficients_S = segment.lineCoefficients(s1, s2);
  Vector2d coefficients_T = segment.lineCoefficients(t1, t2);

  //Costruiamo i segmenti
  Segment segment_S = Segment(s1, s2, coefficients_S);
  Segment segment_T = Segment(t1, t2, coefficients_T);

  //Confrontiamo gli angoli
  double expectedAngleDegree = 90.0;
  double angleDegree = getAngle(segment_S, segment_T);
  EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}

TEST(TestSegment, TestGetAngleWithRandomSegments)
{
  Segment segment;

  //Costruiamo due coppie di punti, per poter costruire i due segmenti
  Point s1 = Point(1.0, 1.0);  //Primo punto, primo segmento
  Point s2 = Point(4.0, 1.0);  //Secondo punto, primo segmento
  Point t1 = Point(1.0, 2.0);  //Primo punto, secondo segmento
  Point t2 = Point(3.0, 4.0);  //Secondo punto, secondo segmento


  //Calcoliamo i coefficienti delle rette associate ai due segmenti
  Vector2d coefficients_S = segment.lineCoefficients(s1, s2);
  Vector2d coefficients_T = segment.lineCoefficients(t1, t2);

  //Costruiamo i segmenti
  Segment segment_S = Segment(s1, s2, coefficients_S);
  Segment segment_T = Segment(t1, t2, coefficients_T);

  //Confrontiamo gli angoli
  double expectedAngleDegree = 45.0;
  double angleDegree = getAngle(segment_S, segment_T);
  EXPECT_DOUBLE_EQ(angleDegree, expectedAngleDegree);
}


///TEST TRIANGLE
TEST(TestTriangle, Test_antiClockWiseOrder)  {

  // Definisco tre punti non ordinati in senso antiorario
    Point p1 = Point(0.0, 0.0);
    Point p2 = Point(1.0, 1.0);
    Point p3 = Point(2.0, 0.0);

  // Creo Triangle
  Triangle triangle = Triangle(p1, p2, p3);

  // Chiamo funzione anticlockwise
  antiClockWiseOrder(triangle);

  // Verifico che i punti siano in senso antiorario
   Point expectedP1 = Point(0.0, 0.0);  //punti giusti in senso antiorario
   Point expectedP2 = Point(2.0, 0.0);
   Point expectedP3 = Point(1.0, 1.0);
   EXPECT_EQ(areSamePoint(triangle.p1, expectedP1), true);
   EXPECT_EQ(areSamePoint(triangle.p2, expectedP2), true);
   EXPECT_EQ(areSamePoint(triangle.p3, expectedP3), true);
}

//CASO 4 ESTREMI DISTINTI
TEST(TestTriangle, TestFirstTriangleWith4DistinctPoints)
{
  Delaunay delaunay;

  // Creo  mesh, un vettore di indici vuoto, un vettore vuoto di triangoli e un ConvexHull vuoto
  vector<Triangle> triangles;
  vector<Segment> convexHull;
  vector<Point> hullPoints;
  vector<int> hullTrianglesIndices;

  // Aggiungo i 6 punti
  vector<Point> points = {Point(1.0, 3.0), Point(4.0, 4.0), Point(0.0, 5.0), Point(5.0, 6.0),  Point(6.0, 2.0),  Point(2.0, 1.0)};

  // Ordino i punti per x crescente e y crescente secondo l'algoritmo Mergesort O(nlogn)
  MergeSorter toBeSortedX;
  MergeSorter toBeSortedY;
  int sx = 0;
  int dx = points.size()-1;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.MergeSort(points, sx, dx, sortingByX);   //(vettore, sx, dx, sortingByX)
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.MergeSort(points, sx, dx, sortingByX);


  // Chiamo funzione firstTriangle
  delaunay.firstTriangle(sortedX, sortedY, triangles, convexHull, hullPoints, hullTrianglesIndices);

  Triangle bestTriangle = triangles[0];
  double bestArea = bestTriangle.Area(bestTriangle.p1, bestTriangle.p2, bestTriangle.p3);

  EXPECT_EQ(bestArea, 11.0);  //triangolo area > è il triangolo 1
}

//CASO 3 ESTREMI DISTINTI
TEST(TestTriangle, TestFirstTriangleWith3DistinctPoints)
{
  Delaunay delaunay;

  // Creo  mesh, un vettore di indici vuoto, un vettore vuoto di triangoli e un ConvexHull vuoto
  vector<Triangle> triangles;
  vector<Segment> convexHull;
  vector<Point> hullPoints;
  vector<int> hullTrianglesIndices;

  // Aggiungo i 6 punti
  vector<Point> points = {Point(1.0, 1.0), Point(2.0, 5.0), Point(3.0, 2.0), Point(4.0, 3.0),  Point(6.0, 1.0),  Point(6.0, 4.0)};
  // Ordino i punti per x crescente e y crescente secondo l'algoritmo Mergesort O(nlogn)
  MergeSorter toBeSortedX;
  MergeSorter toBeSortedY;
  int sx = 0;
  int dx = points.size()-1;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.MergeSort(points, sx, dx, sortingByX);   //(vettore, sx, dx, sortingByX)
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.MergeSort(points, sx, dx, sortingByX);


  // Chiamo funzione firstTriangle
  delaunay.firstTriangle(sortedX, sortedY, triangles, convexHull, hullPoints, hullTrianglesIndices);

  Triangle bestTriangle = triangles[0];
  double bestArea = bestTriangle.Area(bestTriangle.p1, bestTriangle.p2, bestTriangle.p3);

  EXPECT_EQ(bestArea, 8.5);
}

//CASO 2 ESTREMI DISTINTI
TEST(TestTriangle, TestFirstTriangleWith2DistinctPoints)
{
  Delaunay delaunay;

  // Creo  mesh, un vettore di indici vuoto, un vettore vuoto di triangoli e un ConvexHull vuoto
  vector<Triangle> triangles;
  vector<Segment> convexHull;
  vector<Point> hullPoints;
  vector<int> hullTrianglesIndices;

  // Aggiungo i 6 punti
  vector<Point> points = {Point(1, 1), Point(2, 5), Point(3, 2), Point(4, 3),  Point(6, 1),  Point(6, 5) };

  // Ordino i punti per x crescente e y crescente secondo l'algoritmo Mergesort O(nlogn)
  MergeSorter toBeSortedX;
  MergeSorter toBeSortedY;
  int sx = 0;
  int dx = points.size()-1;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.MergeSort(points, sx, dx, sortingByX);   //(vettore, sx, dx, sortingByX)
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.MergeSort(points, sx, dx, sortingByX);


  // Chiamo funzione firstTriangle
  delaunay.firstTriangle(sortedX, sortedY, triangles, convexHull, hullPoints, hullTrianglesIndices);

  Triangle bestTriangle = triangles[0];
  double bestArea = bestTriangle.Area(bestTriangle.p1, bestTriangle.p2, bestTriangle.p3);

  EXPECT_EQ(bestArea, 10.0);
}

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
TEST(TestDelaunay, Test_flip3triangles) {

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

    vector<int> expectedneighbors1 = {2,-1,1};
    vector<int> expectedneighbors2 = {3,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedneighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedneighbors2[i]);
    }
}

//verifichiamo altro caso: adjacentSide1 = 2, adjacentSide2 = 1
TEST(TestDelaunay, Test_flip_triangles) {

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

    vector<int> expectedneighbors1 = {-1,3,1};
    vector<int> expectedneighbors2 = {-1,2,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedneighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedneighbors2[i]);
    }
}

//verifichiamo l'aggiornamento di hulltriangleindices
TEST(TestDelaunay, Test_flip_hull) {

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

    vector<int> expectedneighbors1 = {-1,3,1};
    vector<int> expectedneighbors2 = {2,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(triangles[0].neighbors[i], expectedneighbors1[i]);
        EXPECT_EQ(triangles[1].neighbors[i], expectedneighbors2[i]);
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

    vector<int> expectedneighbors1 = {-1,3,2};
    vector<int> expectedneighbors2 = {4,-1,2};
    vector<int> expectedneighbors3 = {-1,0,1};
    vector<int> expectedneighbors4 = {-1,4,0};
    vector<int> expectedneighbors5 = {-1,1,3};



    for (int i=0; i<3; i++) {
        ASSERT_EQ(triangles[0].neighbors[i], expectedneighbors1[i]);
        ASSERT_EQ(triangles[1].neighbors[i], expectedneighbors2[i]);
        ASSERT_EQ(triangles[2].neighbors[i], expectedneighbors3[i]);
        ASSERT_EQ(triangles[3].neighbors[i], expectedneighbors4[i]);
        ASSERT_EQ(triangles[4].neighbors[i], expectedneighbors5[i]);
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

      vector<int> expectedneighbors1 = {-1, 1, 2};
      vector<int> expectedneighbors2 = {-1, -1, 0};
      vector<int> expectedneighbors3 = {-1, 0, -1};

      for (int i=0; i<3; i++) {
          EXPECT_EQ(triangles[2].neighbors[i], expectedneighbors3[i]);
          EXPECT_EQ(triangles[0].neighbors[i], expectedneighbors1[i]);
          EXPECT_EQ(triangles[1].neighbors[i], expectedneighbors2[i]);
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

#endif // __TEST_EMPTY_H
