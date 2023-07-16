#ifndef __POINT_H
#define __POINT_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Eigen;

namespace ProjectLibrary {

    class Point{

        private:
            int id;
            double x;
            double y;
            
        public:
            //costruttori
            Point() = default;                      //def
            Point(double x, double y);            //x,y
            Point(int id, double x, double y);   //id,x,y

            //getters & setters
            double get_x();
            double get_y();
            int get_id();
            void set_x(double x);
            void set_y(double y);
            void set_id(int id);

            //metodi 
            Point getOrder(Point p, bool sortingByX);
            double getDistance(Point p);
            bool areSamePoint(Point p);

    };
}

#endif // __POINT_H
