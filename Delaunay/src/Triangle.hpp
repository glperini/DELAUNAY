#ifndef __TRIANGLE_H
#define __TRIANGLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

#include "Point.hpp"

using namespace std;
using namespace Eigen;

namespace ProjectLibrary {

    class Triangle{

        private:
            //Point p1;
            //Point p2;
            //Point p3;
            vector<Point> points;
            vector<int> neighbors = {-1, -1, -1};   /*contiene per ogni vicino l'indice del triangolo adiacente e il lato su cui fa adiacenza (lato, triangolo)
                                                                   in prima posizione -> identificativo del lato (lato 1 tra vertice 0-1, lato 2 tra vertice 1-2, lato 3 tra vertice 2-3)
                                                                   in seconda posizione -> INDICE del triangolo vicino
                                                                   di default non ho triangoli vicini per il triangolo grande, quindi avranno tutti e 3 indici pari a -1 */
        public:
            //costruttori
            Triangle() = default; //def
            Triangle(vector<Point> points); //vector<Point>

            //getters & setters
            vector<Point> get_points();
            Point get_point(int index);
            vector<int> get_neighbors();
            int get_neighbor(int index);
            void set_points(vector<Point> points);
            void set_point(Point point, int index);
            void set_neighbors(vector<int> neighbors);
            void set_neighbor(int neighbor, int index);

            //metodi 
            double Area();
            void antiClockWiseOrder();
    };
}

#endif // __TRIANGLE_H
