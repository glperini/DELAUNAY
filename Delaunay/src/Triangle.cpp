#include "Triangle.hpp"

namespace ProjectLibrary {

    //Costruttore p1,p2,p3
    Triangle::Triangle(vector<Point> points):
        points(points)
    {}

    //getters & setters
    vector<Point> Triangle::get_points(){
        return points;
    }

    Point Triangle::get_point(int index){
        return points[index];
    }

    vector<int> Triangle::get_neighbors(){
        return neighbors;
    }
    
    int Triangle::get_neighbor(int index){
        return neighbors[index];
    }
    
    void Triangle::set_points(vector<Point> points){
        this -> points = points;
    }

    void Triangle::set_point(Point point, int index){
        this -> points[index] = point;
    }

    void Triangle::set_neighbors(vector<int> neighbors){
        this -> neighbors = neighbors;
    }

    void Triangle::set_neighbor(int neighbor, int index) {
        this -> neighbors[index] = neighbor;
    }

    //metodi 
    double Triangle::Area(){
        double Area = 0;
        Area =  points[0].get_x() * points[1].get_y() - points[1].get_x() * points[0].get_y() +
                points[1].get_x() * points[2].get_y() - points[2].get_x() * points[1].get_y() +
                points[2].get_x() * points[0].get_y() - points[0].get_x() * points[2].get_y();
        Area = 0.5 * abs(Area);
        return Area;
    }

    void Triangle::antiClockWiseOrder(){
        double crossProduct = (points[1].get_x() - points[0].get_x()) * (points[2].get_y() - points[0].get_y()) - 
                              (points[2].get_x() - points[0].get_x()) * (points[1].get_y() - points[0].get_y());

        //se crossProduct>0 significa che i punti sono ordinati in senso antiorario
        //se crossProduct<0 significa che i punti sono ordinati in senso orario
        if (crossProduct<0){                                              
            Point temporary = points[1];
            points[1]= points[2];       //scambiamo p2 e p3 così ora sono ordinati in senso antiorario
            points[2] = temporary;
        }
    }
}
