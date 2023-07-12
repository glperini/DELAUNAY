#ifndef __EMPTY_H
#define __EMPTY_H

#include <iostream>
#include <fstream>
#include <vector>
//#include "list"
#include "Eigen/Eigen"
//#include "map"
#include <cmath>
//#include <set>
#include <algorithm>

using namespace std;
using namespace Eigen;

namespace ProjectLibrary
{

//DICHIARAZIONE CLASSE Point
class Point{
    public:
        int id;
        double x;
        double y;
        Point() = default;
        Point(double& x, double& y);
        Point(int& id, double& x, double& y);

        //per test
        Point(int id, double x, double y);
        Point(double x, double y);
};

//DICHIARAZIONE CLASSE Segment
class Segment{
   public:
        Point p1;
        Point p2;
        Vector2d coefficients;
        Segment() = default;
        Segment(Point& p1, Point& p2, Vector2d& coefficients);
        //Segment(Point& p1, Point& p2);

        // per test
        Segment(Point p1, Point p2);

        Vector2d lineCoefficients(Point& p1, Point& p2);
    };


//DICHIARAZIONE STRUTTURA TriangularMesh
struct TriangularMesh
    {
    unsigned int NumberOfPoints = 0; ///< number of Cell0D
    vector<Point> points = {}; ///< Cell0D coordinates AND < Cell0D id, (size 2 + size 1) x NumberCell0D (x,y)
    const string fileName;
    const string outputFilePath;
    vector<Segment> segments;

    bool importPoints(TriangularMesh& mesh, const string& fileName);
    bool exportResult(const string& outputFilePath, vector<Segment>& segments);
    };

/*
struct MergeSorter
    {
        int sx = 0;
        int cx = 0;
        int dx = 0;
        bool sortingByX = false;
        vector<Point> v;

        vector<Point> MergeSort(vector<Point>& v,
                                int& sx,
                                int& dx,
                                bool& sortingByX);

        void Merge(vector<Point>& v,
                   int& sx,
                   int& cx,
                   int& dx,
                   bool& sortingByX);       //se cosi non puo funzionare allora inserire  bool& sortingByX come attributo

        Point getOrder(Point& point1, Point& point2, bool& sortingByX);
    };
*/

struct HeapSorter
    {
    int n = 0;
    int i = 0;
    bool sortingByX = false;
    vector<Point> v;

    void heapify(vector<Point>& v, int n, int i, bool sortingByX);
    vector<Point> heapSort(vector<Point>& v, bool sortingByX);
    Point getOrder(Point& point1, Point& point2, bool& sortingByX);
    };

//DICHIARAZIONE CLASSE Triangle
class Triangle
    /*class Triangle:
      PUBLIC
      attributes   =    vector<int> vertices[3], vector<int> neighbors[3]
      constructor  =    default, Triangle(vector<int>& vertices[3], TriangularMesh& mesh)
      methods      =    double Area(vector<int>& vertices[3], TriangularMesh& mesh) */
    {
    public:
        Point p1;
        Point p2;
        Point p3;

        vector<int> neighbors ={-1, -1, -1};   //contiene per ogni vicino l'indice del triangolo adiacente e il lato su cui fa adiacenza (lato, triangolo)
                                                                   //in prima posizione -> identificativo del lato (lato 1 tra vertice 0-1, lato 2 tra vertice 1-2, lato 3 tra vertice 2-3)
                                                                   //in seconda posizione -> INDICE del triangolo vicino
                                                                   //di default non ho triangoli vicini per il triangolo grande, quindi avranno tutti e 3 indici pari a -1
        Triangle() = default;
        Triangle(Point& p1, Point& p2, Point& p3);

        double Area(Point& p1, Point& p2, Point& p3);
        //void antiClockWiseOrder(Point& p1, Point& p2, Point& p3);

        //per test
        //Triangle(Point p1, Point p2, Point p3);
    };


//DICHIARAZIONE CLASSE Delayunay
class Delaunay
    {
    /*class Delaunay:
      PUBLIC
      attributes   =    vector<vector> coordinates, vector<Triangle> triangles
      constructor  =    DEFAULT, Delaunay(const vector<vector>& coordinates),
      methods      =    vector<Triangle> getTriangulation(),
                        bool importPoints(TriangularMesh& mesh, string& fileName, vector& leftOutIndices),
                        void void firstTriangle(TriangularMesh& mesh, vector& indices, vector<Triangle>& triangles),
                        void insertPoint(vector& indices, vector& leftOutIndices),

                        int findTriangleContainingPoint(const vector& point),
                        void verifyDelaunayCondition(const Triangle& triangle, int pointIndex),
                        bool isPointInsideCircumcircle(const vector& point, const Triangle& triangle),
                        void flipTriangles*/
    public:
        vector<Point> points;           //vettore di punti
        vector<Triangle> triangles;     //vettore che contiene i triangoli aggiunti alla triangolazione
        vector<Segment> convexHull;     //guscio esterno
        vector<Point> hullPoints;       //punti del guscio
        vector<int> hullTrianglesIndices; //indici di triangoli del guscio
        int triangleIndex = 0;          //indice di un triangolo che stiamo testando

        //Point point;                    //punto che stiamo considerando

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

        /*void onSegment(Point& point,                                    //il nuovo punto si trova su un lato già esistente
                       int& triangleIndex,
                       int& sideWithPoint,
                       vector<Triangle>& triangles,
                       vector<Segment>& convexHull,
                       vector<int>& hullTrianglesIndices);*/

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


//DICHIARAZIONE FUNZIONI AUSILIARIE
///DEF funzione
inline int binarySearch(vector<Point> points, Point target){
    int sx = 0;
    int dx = points.size();
    bool found = false;
    int index = -1;
    while ((sx<=dx) && (!found)){
        int cx = (sx+dx)/2;
        if (points[cx].x>target.x)
            dx = cx - 1;
        else if (points[cx].x<target.x)
            sx = cx + 1;
        else {
            if (points[cx].y==target.y){
                index = cx;
                found = true;
            }
            else if (points[cx].y>target.y){
                cx--;
                while (!found){
                    if (points[cx].y==target.y){
                        index = cx;
                        found = true;
                    }
                    cx--;
                }
            }
            else {
                cx++;
                while (!found){
                    if (points[cx].y==target.y){
                        index = cx;
                        found = true;
                    }
                    cx++;
                }
            }
        }
    }
    return index;
}



///DEF funzione
//Calcola la distanza fra due punti
inline double getDistance(Point p1, Point p2)
{
  return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

///DEF funzione
inline Point findIntersection(Vector2d r1, Vector2d r2){  //qui un vettore è formato da m e q
    Point intersectionPoint;
    if (!isnan(r1[0]) && !isnan(r2[0])){  //se le rette non sono verticali
       if (r1[0]==r2[0]){
           if (r1[1]==r2[1])        //caso di due rette coincidenti
               intersectionPoint = Point(10, r1[0]*10 + r1[1]);     //assegnato valore casuale solo per avere punto valido
           else
               intersectionPoint = Point(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
       }
       else {      // (r1[0]!=r2[0])
           double x = (r2[1]-r1[1])/(r1[0]-r2[0]);     //si trova la x con x = (q2-q1)/(m1-m2)
           intersectionPoint = Point(x, r1[0]*x + r1[1]);     //si sostituisce in una delle due rette per trovare la y
       }
   }
   //se invece una retta è verticale o entrambe sono verticali:
   else if (isnan(r1[0]) && isnan(r2[0]))       //se entrambe le rette sono verticali non c'è intersezione
       intersectionPoint = Point(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
   else if (isnan(r1[0]))     //se la prima retta è verticale
       intersectionPoint = Point(r1[1], r2[0]*r1[1] + r2[1]);
   else if (isnan(r2[0]))     //se la seconda retta è verticale
       intersectionPoint = Point(r2[1], r1[0]*r2[1] + r1[1]);
   return intersectionPoint;
}

///DEF funzione
//Calcola angolo tra due segmenti
inline double getAngle(Segment segment1, Segment segment2){
    double angleDegree;
    //verifico se i due segmenti hanno stesso m (cioè sono paralleli) oppure se sono entrambi verticali
    if (segment1.coefficients[0]==segment2.coefficients[0] || (isnan(segment1.coefficients[0]) && isnan(segment2.coefficients[0])))
        angleDegree = 0.0;
    //se i due segmenti sono perpendicolari (ossia m1*m2 = -1):
    else if (segment1.coefficients[0]*segment2.coefficients[0]==-1)
        angleDegree = 90.0;
    //se i due segmenti non sono nè paralleli nè perpendicolari:
    else {
        double vector1X = segment1.p2.x - segment1.p1.x;
        double vector1Y = segment1.p2.y - segment1.p1.y;
        double vector2X = segment2.p2.x - segment2.p1.x;
        double vector2Y = segment2.p2.y - segment2.p1.y;

        // Calcola i prodotti scalari dei due vettori
        double dotProduct = vector1X * vector2X + vector1Y * vector2Y;

        // Calcola le lunghezze dei due vettori
        double vector1Length = sqrt(vector1X * vector1X + vector1Y * vector1Y);
        double vector2Length = sqrt(vector2X * vector2X + vector2Y * vector2Y);

        // Calcola l'angolo tra i due vettori usando la formula dell'arcocoseno
        double angle = acos(dotProduct / (vector1Length * vector2Length));

        //converto l'angolo in gradi
        angleDegree = angle * 180.0 / M_PI;
    }
    return angleDegree;
}

///DEF funzione di uguaglianza tra vettori
template<typename T>
inline bool areEqual(T v1, T v2) {
    if (isnan(v1[0]) && isnan(v2[0]) && v1[1] == v2[1])
            return true;
    else if (v1[0] == v2[0] && v1[1] == v2[1])
        return true;
    else
        return false;
}

///DEF funzione di uguaglianza tra punti
inline bool areSamePoint(Point p1, Point p2) {
    if (p1.x == p2.x && p1.y == p2.y)
        return true;
    else
        return false;
}

///DEF funzione di uguaglianza tra segmenti
inline bool areSameSegment(Segment s1, Segment s2) {
    if ((areSamePoint(s1.p1, s2.p1) && areSamePoint(s1.p2, s2.p2)) || (areSamePoint(s1.p1, s2.p2) && areSamePoint(s1.p2, s2.p1)))
        return true;
    else
        return false;
}

///DEF funzione
inline void antiClockWiseOrder(Triangle& triangle)
{
    double crossProduct = (triangle.p2.x - triangle.p1.x) * (triangle.p3.y - triangle.p1.y) - (triangle.p3.x - triangle.p1.x) * (triangle.p2.y - triangle.p1.y);
    if (crossProduct>0)                 //se crossProduct>0 significa che i punti sono ordinati in senso antiorario
        return;
    else {                               //se crossProduct<0 significa che i punti sono ordinati in senso orario
        Point temporary = triangle.p2;
        triangle.p2 = triangle.p3;       //scambaimo p2 e p3 così ora sono ordinati in senso antiorario
        triangle.p3 = temporary;
        return;
    }
}

}

#endif // __EMPTY_H
