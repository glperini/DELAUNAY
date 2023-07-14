#ifndef __OBJECTS_H
#define __OBJECTS_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>


using namespace std;
using namespace Eigen;

namespace LibraryObjects
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

//DICHIARAZIONE CLASSE Triangle
class Triangle
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
    };

}

#endif // __OBJECTS_H
