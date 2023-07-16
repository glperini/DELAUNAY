#ifndef __TEST_DELAUNAY_H
#define __TEST_DELAUNAY_H

#include <gtest/gtest.h>
#include "Delaunay.hpp"
#include "TriangularMesh.hpp"
#include "HeapSort.hpp"

using namespace ProjectLibrary;
using namespace SortingLibrary;

using namespace testing;

///TEST DELAUNAY

//CASO 4 ESTREMI DISTINTI
TEST(TestDelaunay, Test_firstTriangle_4DistinctPoints)
{
    Delaunay delaunay;
    TriangularMesh mesh;

    // Aggiungo i 6 punti
    mesh.set_points({Point({1.0, 3.0}), Point({4.0, 4.0}), Point({0.0, 5.0}), Point({5.0, 6.0}),  Point({6.0, 2.0}),  Point({2.0, 1.0})});

    // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
    bool sortingByX = true;
    vector<Point> sortedX = heapSort(mesh.get_points(), sortingByX);
    sortingByX = !sortingByX;
    vector<Point> sortedY = heapSort(mesh.get_points(), sortingByX);


    // Chiamo funzione firstTriangle
    delaunay.firstTriangle(sortedX, sortedY);

    double bestArea = delaunay.get_triangle(0).Area();

    EXPECT_EQ(bestArea, 11.0);  //triangolo area > è il triangolo 1
}

//CASO 3 ESTREMI DISTINTI
TEST(TestDelaunay, Test_firstTriangle_3DistinctPoints)
{
    Delaunay delaunay;
    TriangularMesh mesh;

    // Aggiungo i 6 punti
    mesh.set_points({Point({1.0, 1.0}), Point({2.0, 5.0}), Point({3.0, 2.0}), Point({4.0, 3.0}),  Point({6.0, 1.0}),  Point({6.0, 4.0})});

    // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
    bool sortingByX = true;
    vector<Point> sortedX = heapSort(mesh.get_points(), sortingByX);
    sortingByX = !sortingByX;
    vector<Point> sortedY = heapSort(mesh.get_points(), sortingByX);

    // Chiamo funzione firstTriangle
    delaunay.firstTriangle(sortedX, sortedY);

    double bestArea = delaunay.get_triangle(0).Area();

    EXPECT_EQ(bestArea, 8.5);
}

//CASO 2 ESTREMI DISTINTI
TEST(TestDelaunay, Test_firstTriangle_2DistinctPoints)
{
    Delaunay delaunay;
    TriangularMesh mesh;

    // Aggiungo i 6 punti
    mesh.set_points({Point({1, 1}), Point({2, 5}), Point({3, 2}), Point({4, 3}),  Point({6, 1}),  Point({6, 5})});

    // Ordino i punti per x crescente e y crescente secondo l'algoritmo HeapSort O(nlogn)
    bool sortingByX = true;
    vector<Point> sortedX = heapSort(mesh.get_points(), sortingByX);
    sortingByX = !sortingByX;
    vector<Point> sortedY = heapSort(mesh.get_points(), sortingByX);

    // Chiamo funzione firstTriangle
    delaunay.firstTriangle(sortedX, sortedY);

    double bestArea = delaunay.get_triangle(0).Area();

    EXPECT_EQ(bestArea, 10.0);
}

TEST(TestDelaunay, Test_isPointInsideCircumcircle_In) {
    Delaunay delaunay;

    // Creo  triangoli di prova
    Point p1 = Point({0.0, 0.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({0.0, 3.0});
    Point p4 = Point({3.0, 3.0});
    Point p5 = Point({5.0, 0.0});

    delaunay.set_triangles({Triangle({p1, p2, p3}), Triangle({p2, p4, p3}), Triangle({p2, p5, p4})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideCircle = Point({4.0, 1.0});

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideCircumcircle(pointInsideCircle);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, 2);
}

TEST(TestDelaunay, Test_isPointInsideCircumcircle_Out) {
    Delaunay delaunay;

    // Creo  triangoli di prova
    Point p1 = Point({0.0, 0.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({0.0, 3.0});
    Point p4 = Point({3.0, 3.0});
    Point p5 = Point({5.0, 0.0});

    delaunay.set_triangles({Triangle({p1, p2, p3}), Triangle({p2, p4, p3}), Triangle({p2, p5, p4})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideCircle = Point({10.0, 0.0});

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideCircumcircle(pointInsideCircle);

    // Verifico che la funzione abbia restituito un valore diverso da -1
    EXPECT_EQ(resultIndex, -1);
}


TEST(TestDelaunay, Test_isPointInsideTriangle_In) {
    Delaunay delaunay;

    // Creo  triangoli di prova
    Point p1 = Point({0.0, 0.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({0.0, 3.0});

    delaunay.set_triangles({Triangle({p1, p2, p3})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideTriangle = Point({1.0, 1.0});

    //chiamo la funzione
    int resultIndex = delaunay.isPointInsideTriangle(pointInsideTriangle);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, 0);
}

TEST(TestDelaunay, Test_isPointInsideTriangle_Out) {

    Delaunay delaunay;

    // Creo  triangoli di prova
    Point p1 = Point({0.0, 0.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({0.0, 3.0});

    delaunay.set_triangles({Triangle({p1, p2, p3})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointInsideTriangle = Point({2.0, 2.0});

    //chiamo la funzione
   int resultIndex = delaunay.isPointInsideTriangle(pointInsideTriangle);

    // Verifica che l'indice del triangolo sia 0
    EXPECT_EQ(resultIndex, -1);
}

TEST(TestDelaunay, Test_findTriangleContainingPoint_In) {

    Delaunay delaunay;

    //Creo triangoli di prova
    Point p1 = Point({2.0, 1.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({5.0, 1.5});
    Point p4 = Point({4.0, 2.0});
    Point p5 = Point({5.0, 0.0});

    delaunay.set_triangles({Triangle({p1, p2, p3}), Triangle({p2, p5, p3}), Triangle({p1, p3, p4})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointOutsideTriangle = Point({4.0, 1.5});

    //chiamo la funzione
    int resultIndex = delaunay.findTriangleContainingPoint(pointOutsideTriangle);

    EXPECT_EQ(resultIndex, 2);
}

TEST(TestDelaunay, Test_findTriangleContainingPoint_Out) {

    Delaunay delaunay;

    //Creo triangoli di prova
    Point p1 = Point({2.0, 1.0});
    Point p2 = Point({3.0, 0.0});
    Point p3 = Point({5.0, 1.5});
    Point p4 = Point({4.0, 2.0});
    Point p5 = Point({5.0, 0.0});

    delaunay.set_triangles({Triangle({p1, p2, p3}), Triangle({p2, p5, p3}), Triangle({p1, p3, p4})});

    // Testo un punto che si trova all'interno del circumcerchio del triangolo
    Point pointOutsideTriangle = Point({2.0, 0.0});

    //chiamo la funzione
    int resultIndex = delaunay.findTriangleContainingPoint(pointOutsideTriangle);

    EXPECT_EQ(resultIndex, -1);
}

//verifichiamo caso: adjacentSide1 = 1, adjacentSide2 = 2
TEST(TestDelaunay, Test_flipTriangles_Side12Neighbors) {

    Delaunay delaunay;

    Point p1 = Point({1.0, 2.0});
    Point p2 = Point({2.0, 0.0});
    Point p3 = Point({2.0, 4.0});
    Point p4 = Point({3.5, 2.0});
    Point p5 = Point({1.0, 0.0});
    Point p6 = Point({3.5, 4.5});
    Triangle triangle1 = Triangle({p1,p2,p3});
    Triangle triangle2 = Triangle({p2,p4,p3});
    Triangle triangle3 = Triangle({p1,p5,p2});
    Triangle triangle4 = Triangle({p4,p6,p3});
    delaunay.set_triangles({triangle1, triangle2, triangle3, triangle4});

    int adjacentSide1 = 1;
    int adjacentSide2 = 2;
    delaunay.set_hullTrianglesIndices({0, 1, 2, 3});
    int triangle1Index = 0;
    int triangle2Index = 1;

    Triangle t0 = delaunay.get_triangle(0);
    Triangle t1 = delaunay.get_triangle(1);
    Triangle t2 = delaunay.get_triangle(2);
    Triangle t3 = delaunay.get_triangle(3);

    t0.set_neighbors({2,1,-1});
    t1.set_neighbors({-1,3,0});
    t2.set_neighbors({-1,-1,0});
    t3.set_neighbors({-1,-1,1});

    delaunay.set_triangle(t0,0);
    delaunay.set_triangle(t1,1);
    delaunay.set_triangle(t2,2);
    delaunay.set_triangle(t3,3);

    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2);

    vector<int> expectedNeighbors1 = {2,-1,1};
    vector<int> expectedNeighbors2 = {3,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(delaunay.get_triangle(0).get_neighbor(i), expectedNeighbors1[i]);
        EXPECT_EQ(delaunay.get_triangle(1).get_neighbor(i), expectedNeighbors2[i]);
    }
}

//verifichiamo altro caso: adjacentSide1 = 2, adjacentSide2 = 1
TEST(TestDelaunay, Test_flipTriangles_Side21Neighbors) {

    Delaunay delaunay;

    Point p1 = Point({1.0, 2.0});
    Point p2 = Point({2.0, 0.0});
    Point p3 = Point({2.0, 4.0});
    Point p4 = Point({3.5, 2.0});
    Point p5 = Point({1.0, 0.0});
    Point p6 = Point({3.5, 4.5});

    Triangle triangle1 = Triangle({p1,p2,p4});
    Triangle triangle2 = Triangle({p3,p1,p4});
    Triangle triangle3 = Triangle({p1,p5,p2});
    Triangle triangle4 = Triangle({p4,p6,p3});
    delaunay.set_triangles({triangle1, triangle2, triangle3, triangle4});

    int adjacentSide1 = 2;
    int adjacentSide2 = 1;
    delaunay.set_hullTrianglesIndices({0, 1, 2, 3});
    int triangle1Index = 0;
    int triangle2Index = 1;

    Triangle t0 = delaunay.get_triangle(0);
    Triangle t1 = delaunay.get_triangle(1);
    Triangle t2 = delaunay.get_triangle(2);
    Triangle t3 = delaunay.get_triangle(3);

    t0.set_neighbors({2,-1,1});
    t1.set_neighbors({-1,0,3});
    t2.set_neighbors({-1,-1,0});
    t3.set_neighbors({-1,-1,1});

    delaunay.set_triangle(t0,0);
    delaunay.set_triangle(t1,1);
    delaunay.set_triangle(t2,2);
    delaunay.set_triangle(t3,3);

    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2);

    vector<int> expectedNeighbors1 = {-1,3,1};
    vector<int> expectedNeighbors2 = {-1,2,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(delaunay.get_triangle(0).get_neighbor(i), expectedNeighbors1[i]);
        EXPECT_EQ(delaunay.get_triangle(1).get_neighbor(i), expectedNeighbors2[i]);
    }
}

//verifichiamo l'aggiornamento di hullTriangleIndices
TEST(TestDelaunay, Test_flipTriangles_hullTriangleIndices) {

    Delaunay delaunay;

    Point p1 = Point({1.0, 2.0});
    Point p2 = Point({2.0, 0.0});
    Point p3 = Point({2.0, 4.0});
    Point p4 = Point({3.5, 2.0});
    Point p5 = Point({1.0, 0.0});
    Point p6 = Point({3.5, 4.5});
    Point p7 = Point({3.0, 0.0});

    Triangle triangle1 = Triangle({p1,p2,p3});
    Triangle triangle2 = Triangle({p2,p4,p3});
    Triangle triangle3 = Triangle({p1,p5,p2});
    Triangle triangle4 = Triangle({p4,p6,p3});
    Triangle triangle5 = Triangle({p2,p7,p4});
    delaunay.set_triangles({triangle1, triangle2, triangle3, triangle4, triangle5});

    int adjacentSide1 = 1;
    int adjacentSide2 = 2;
    delaunay.set_hullTrianglesIndices({0, 2, 3, 4});
    int triangle1Index = 0;
    int triangle2Index = 1;

    Triangle t0 = delaunay.get_triangle(0);
    Triangle t1 = delaunay.get_triangle(1);
    Triangle t2 = delaunay.get_triangle(2);
    Triangle t3 = delaunay.get_triangle(3);
    Triangle t4 = delaunay.get_triangle(4);

    t0.set_neighbors({2,1,-1});
    t1.set_neighbors({4,3,0});
    t2.set_neighbors({-1,-1,0});
    t3.set_neighbors({-1,-1,1});
    t4.set_neighbors({-1,-1,1});

       //idea per risolvere l'avviso
    delaunay.set_triangle(t0,0);
    delaunay.set_triangle(t1,1);
    delaunay.set_triangle(t2,2);
    delaunay.set_triangle(t3,3);


    delaunay.flipTriangles(triangle1Index, triangle2Index, adjacentSide1, adjacentSide2);

    vector<int> expectedHullIndices = {2,3,4,1};
    for (int i=0;i<4;i++)
        EXPECT_EQ(delaunay.get_hullTrianglesIndex(i), expectedHullIndices[i]);
}

//
TEST(TestDelaunay, Test_verifyDelaunayCondition) {

    Delaunay delaunay;

    Point p1 = Point({1.0, 2.0});
    Point p2 = Point({2.0, 0.0});
    Point p3 = Point({2.0, 4.0});
    Point p4 = Point({3.5, 2.0});
    Point p5 = Point({3.5, 4.5});
    Point p6 = Point({3.0, 0.0});

    delaunay.set_triangleIndex(1);

    Triangle triangle1 = Triangle({p1,p2,p3});
    Triangle triangle2 = Triangle({p2,p4,p3});
    Triangle triangle3 = Triangle({p4,p5,p3});
    Triangle triangle4 = Triangle({p2,p6,p4});
    delaunay.set_triangles({triangle1, triangle2, triangle3, triangle4});

    delaunay.set_hullTrianglesIndices({0,2,3});
    Triangle t0 = delaunay.get_triangle(0);
    Triangle t1 = delaunay.get_triangle(1);
    Triangle t2 = delaunay.get_triangle(2);
    Triangle t3 = delaunay.get_triangle(3);

    t0.set_neighbors({-1,1,-1});
    t1.set_neighbors({3,2,0});
    t2.set_neighbors({-1,-1,1});
    t3.set_neighbors({-1,-1,1});

    delaunay.set_triangle(t0,0);
    delaunay.set_triangle(t1,1);
    delaunay.set_triangle(t2,2);
    delaunay.set_triangle(t3,3);

    delaunay.verifyDelaunayCondition(delaunay.get_triangle(1));

    vector<int> expectedNeighbors1 = {-1,3,1};
    vector<int> expectedNeighbors2 = {2,-1,0};

    for (int i=0; i<3; i++) {
        EXPECT_EQ(delaunay.get_triangle(0).get_neighbor(i), expectedNeighbors1[i]);
        EXPECT_EQ(delaunay.get_triangle(1).get_neighbor(i), expectedNeighbors2[i]);
    }

    vector<int> expectedHullTriangleIndices = {0,2,3,1};
    for (int i=0; i<4; i++)
        EXPECT_EQ(delaunay.get_hullTrianglesIndex(i), expectedHullTriangleIndices[i]);
}

TEST(TestDelaunay, Test_inTriangle) {
    Delaunay delaunay;

    Point p1 = Point({3.0, 3.0});
    Point p2 = Point({5.0, 7.0});
    Point p3 = Point({11.0, 3.0});
    Point p4 = Point({6.0, 2.0});
    Point p5 = Point({3.0, 6.0});
    Point InsidePoint = Point({6.0, 5.0});

    delaunay.set_triangleIndex(1);

    Triangle triangle1 = Triangle({p1,p4,p2});
    Triangle triangle2 = Triangle({p2,p4,p3});
    Triangle triangle3 = Triangle({p1,p2,p5});
    delaunay.set_triangles({triangle1, triangle2, triangle3});

    delaunay.set_hullTrianglesIndices({0, 1, 2});

    Triangle t0 = delaunay.get_triangle(0);
    Triangle t1 = delaunay.get_triangle(1);
    Triangle t2 = delaunay.get_triangle(2);

    t0.set_neighbors({-1,1,2});
    t1.set_neighbors({0,-1,-1});
    t2.set_neighbors({0,-1,-1});

    delaunay.set_triangle(t0,0);
    delaunay.set_triangle(t1,1);
    delaunay.set_triangle(t2,2);

    delaunay.inTriangle(InsidePoint);

    vector<int> expectedNeighbors1 = {-1,3,2};
    vector<int> expectedNeighbors2 = {4,-1,2};
    vector<int> expectedNeighbors3 = {-1,0,1};
    vector<int> expectedNeighbors4 = {-1,4,0};
    vector<int> expectedNeighbors5 = {-1,1,3};


    for (int i=0; i<3; i++) {
        ASSERT_EQ(delaunay.get_triangle(0).get_neighbor(i), expectedNeighbors1[i]);
        ASSERT_EQ(delaunay.get_triangle(1).get_neighbor(i), expectedNeighbors2[i]);
        ASSERT_EQ(delaunay.get_triangle(2).get_neighbor(i), expectedNeighbors3[i]);
        ASSERT_EQ(delaunay.get_triangle(3).get_neighbor(i), expectedNeighbors4[i]);
        ASSERT_EQ(delaunay.get_triangle(4).get_neighbor(i), expectedNeighbors5[i]);
    }

}

TEST(TestDelaunay, Test_outTriangle) {

      Delaunay delaunay;

      Point outsidePoint = Point({-2.0, 2.0});

      Point p1 = Point({0.0, 0.0});
      Point p2 = Point({3.0, 0.0});
      Point p3 = Point({1.0, 2.0});
      Point p4 = Point({4.0, 2.0});

      delaunay.set_triangles({Triangle({p1, p2, p3}), Triangle({p2, p4, p3})});


      Triangle t0 = delaunay.get_triangle(0);
      Triangle t1 = delaunay.get_triangle(1);

      t0.set_neighbors({-1,1,-1});
      t1.set_neighbors({-1,-1,0});

      delaunay.set_triangle(t0,0);
      delaunay.set_triangle(t1,1);

      Segment temp;
      delaunay.set_convexHull({Segment({p1, p2}), Segment({p2, p4}), Segment({p4, p3}), Segment({p3, p1})});
      temp = delaunay.get_convexHull_s(0);  //{0.0, 0.0}
      temp.lineCoefficients();
      delaunay.set_convexHull_s(temp,0);

      temp = delaunay.get_convexHull_s(1);  //{2.0, -6.0}
      temp.lineCoefficients();
      delaunay.set_convexHull_s(temp,1);

      temp = delaunay.get_convexHull_s(2);  //{0.0, 2.0}
      temp.lineCoefficients();
      delaunay.set_convexHull_s(temp,2);

      temp = delaunay.get_convexHull_s(3);  //{2.0, 0.0}
      temp.lineCoefficients();
      delaunay.set_convexHull_s(temp,3);

      delaunay.set_hullPoints({p1, p2, p3, p4});

      delaunay.set_hullTrianglesIndices({0, 1});

      delaunay.outTriangle(outsidePoint);

      //verifico che il triangolo creato sia corretto
      EXPECT_EQ(delaunay.get_triangle(2).get_point(0).areSamePoint(outsidePoint), true);
      EXPECT_EQ(delaunay.get_triangle(2).get_point(1).areSamePoint(p1), true);
      EXPECT_EQ(delaunay.get_triangle(2).get_point(2).areSamePoint(p3), true);

      vector<int> expectedNeighbors1 = {-1, 1, 2};
      vector<int> expectedNeighbors2 = {-1, -1, 0};
      vector<int> expectedNeighbors3 = {-1, 0, -1};

      for (int i=0; i<3; i++) {
          EXPECT_EQ(delaunay.get_triangle(2).get_neighbor(i), expectedNeighbors3[i]);
          EXPECT_EQ(delaunay.get_triangle(0).get_neighbor(i), expectedNeighbors1[i]);
          EXPECT_EQ(delaunay.get_triangle(1).get_neighbor(i), expectedNeighbors2[i]);
          }

      //verifico i punti nel guscio: expectedHullPoints = {p1, p2, p3, p4, outsidePoint};
      EXPECT_EQ(delaunay.get_hullPoint(4).areSamePoint(outsidePoint), true);

      //verifico che abbia aggiunto il nuovo triangolo nel guscio: expectedHullTrianglesIndices = {0, 1, 2}
      EXPECT_EQ(delaunay.get_hullTrianglesIndex(2), 2);
}

#endif // __TEST_DELAUNAY_H
