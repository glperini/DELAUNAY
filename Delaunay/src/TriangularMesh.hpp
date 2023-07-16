#ifndef __TRIANGULARMESH_H
#define __TRIANGULARMESH_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

#include "Triangle.hpp"
#include "Segment.hpp"

using namespace std;
using namespace Eigen;

namespace ProjectLibrary {

    class TriangularMesh{

        private:
            unsigned int numberOfPoints = 0; ///< number of Cell0D
            vector<Point> points; ///< Cell0D coordinates AND < Cell0D id, (size 2 + size 1) x NumberCell0D (x,y)
            string fileName;
            string outputFilePath;
            vector<Triangle> triangles {};
            vector<Segment> segments = {};

        public:
            //costruttore
            TriangularMesh() = default; //def

            //getters & setters
            unsigned int get_numberOfPoints();
            vector<Point> get_points();
            const string get_fileName();
            const string get_outputFilePath();
            vector<Triangle> get_triangles();
            vector<Segment> get_segments();
            Segment get_segment(int index);
            void set_numberOfPoints(unsigned int numberOfPoints);
            void set_points(vector<Point> points);
            void set_fileName(const string fileName);
            void set_outputFilePath(string outputFilePath);
            void set_triangles(vector<Triangle> triangles);
            void set_segments(vector<Segment> segments);
            void set_segment(Segment segment, int index);

            //metodi 
            bool importPoints();
            void drawSegments();
            bool exportResult();
    
    };
}

#endif // __TRIANGULARMESH_H
