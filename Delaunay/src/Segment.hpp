#ifndef __SEGMENT_H
#define __SEGMENT_H

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

    class Segment{

        private:
            vector<Point> points;
            Vector2d coefficients;

        public:
            //costruttori
            Segment() = default;                    //def
            Segment(vector<Point> points);       //points

            //getters & setters
            vector<Point> get_points();
            Point get_point(int index);
            Vector2d get_coefficients();
            double get_coefficient(int index);
            void set_points(vector<Point> points);
            void set_point(Point point, int index);
            void set_coefficients(Vector2d coefficients);
            void set_coefficient(double coefficient, int index);

            //metodi 
            void lineCoefficients();
            double getAngle(Segment segment);
            bool areSameSegment(Segment segment);
    };
}

#endif // __SEGMENT_H
