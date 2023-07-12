#ifndef __DELAUNAY_H
#define __DELAUNAY_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

#include "sorting.hpp"
#include "functions.hpp"

using namespace LibraryObjects;
using namespace LibrarySorting;
using namespace LibraryFunctions;

using namespace std;
using namespace Eigen;

namespace LibraryDelaunay
{

//DICHIARAZIONE CLASSE Delayunay
class Delaunay
    {
    public:
        vector<Point> points;           //vettore di punti
        vector<Triangle> triangles;     //vettore che contiene i triangoli aggiunti alla triangolazione
        vector<Segment> convexHull;     //guscio esterno
        vector<Point> hullPoints;       //punti del guscio
        vector<int> hullTrianglesIndices; //indici di triangoli del guscio
        int triangleIndex = 0;          //indice di un triangolo che stiamo testando

        Delaunay() = default;

        vector<Triangle> getTriangulation();


        void firstTriangle(vector<Point>& sortedX,                       //per creare il primo triangolo
                           vector<Point>& sortedY,
                           vector<Triangle>& triangles,
                           vector<Segment>& convexHull,
                           vector<Point>& hullPoints,
                           vector<int>& hullTrianglesIndices);

        int isPointInsideCircumcircle(Point& point,                     //se punto è interno cerchio circoscritto al triangolo: restituisce indice triangolo se trovato, -1 altrimenti
                                      vector<Triangle>& triangles,
                                      int& triangleIndex);


        int isPointInsideTriangle(Point& point,                    //per verificare se il punto appena aggiunto si trova dentro o fuori alla triangolazione e agire di conseguenza
                                  vector<Triangle>& triangles,
                                  int& triangleIndex);

        int findTriangleContainingPoint(Point& point,              //cerca il triangolo contenente un determinato punto
                                        vector<Triangle>& triangles,
                                        int& triangleIndex);

        void inTriangle(Point& point,                                   //il nuovo punto inserito si trova all'interno della triangolazione
                        vector<Triangle>& triangles,
                        int& triangleIndex,
                        vector<int>& hullTrianglesIndices);

        void outTriangle(Point& outsidePoint,                           //il nuovo punto inserito si trova all'esterno della triangolazione
                         vector<Triangle>& triangles,
                         vector<Segment>& convexHull,
                         vector<Point>& hullPoints,
                         vector<int>& hullTrianglesIndices);

        void verifyDelaunayCondition(Triangle& triangle,                //verifica l'ipotesi di Delaunay
                                     int& triangleIndex,
                                     vector<Triangle>& triangles,
                                     vector<int>& hullTrianglesIndices);

        void flipTriangles(int& triangle1Index,   //svolge operazione di flip
                           int& triangle2Index,
                           int& adjacentSide1,     //rispetto al triangolo 1, il triangolo 2 adiacente che sto considerando si trova su questo lato
                           int& adjacentSide2,     //rispetto al triangolo 2, il triangolo 1 adiacente che sto considerando si trova su questo lato
                           vector<int>& hullTrianglesIndices,
                           vector<Triangle>& triangles);

        vector<Segment> drawSegments(vector<Triangle>& triangles);

    };

}

#endif // __DELAUNAY_H
