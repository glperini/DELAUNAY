#ifndef __TEST_OBJECTS_H
#define __TEST_OBJECTS_H

#include <gtest/gtest.h>
#include "objects.hpp"
#include "sorting.hpp"
#include "functions.hpp"
#include "delaunay.hpp"

using namespace LibraryObjects;
using namespace LibrarySorting;
using namespace LibraryDelaunay;
using namespace LibraryFunctions;

using namespace testing;

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

  // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
  HeapSorter toBeSortedX;
  HeapSorter toBeSortedY;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.heapSort(points, sortingByX);
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.heapSort(points, sortingByX);


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
  // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
  HeapSorter toBeSortedX;
  HeapSorter toBeSortedY;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.heapSort(points, sortingByX);
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.heapSort(points, sortingByX);


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

  // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
  HeapSorter toBeSortedX;
  HeapSorter toBeSortedY;
  bool sortingByX = true;
  vector<Point> sortedX = toBeSortedX.heapSort(points, sortingByX);
  sortingByX = !sortingByX;
  vector<Point> sortedY = toBeSortedY.heapSort(points, sortingByX);


  // Chiamo funzione firstTriangle
  delaunay.firstTriangle(sortedX, sortedY, triangles, convexHull, hullPoints, hullTrianglesIndices);

  Triangle bestTriangle = triangles[0];
  double bestArea = bestTriangle.Area(bestTriangle.p1, bestTriangle.p2, bestTriangle.p3);

  EXPECT_EQ(bestArea, 10.0);
}

#endif // __TEST_OBJECTS_H
