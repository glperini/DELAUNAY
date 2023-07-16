#include "Point.hpp"

namespace ProjectLibrary {

    //Costruttore x,y
    Point::Point(double x, double y):
        x(x),
        y(y)
    {id=0;}

    //Costruttore id,x,y
    Point::Point(int id, double x, double y):
        id(id),
        x(x),
        y(y)
    {}


    //getters & setters
    double Point::get_x(){
        return x;
    }

    double Point::get_y(){
        return y;
    }

    int Point::get_id(){
        return id;
    }

    void Point::set_x(double x){
        this -> x = x;
    }

    void Point::set_y(double y){
        this -> y = y;
    }
    
    void Point::set_id(int id){
        this -> id = id;
    }

    //metodi 
    double Point::getDistance(Point p)
    {
        return sqrt((this->x - p.x)*(this->x - p.x) + (this->y - p.y)*(this->y - p.y));
    }

    bool Point::areSamePoint(Point p){
        if (this->x == p.x && this->y == p.y)
            return true;
        else
            return false;
    }

    Point Point::getOrder(Point p, bool sortingByX){
        Point temporary;
        if (sortingByX == true) {
            if (x<p.get_x()){
                temporary = *this;
                return temporary;
            }
            else if (x>p.get_x()) {
                temporary = p;
                return temporary;
            }
            else { // point1.x == point2.x
                if (y<p.get_y()) {
                    temporary = *this;
                    return temporary;
                }
                else {
                    temporary = p;
                    return temporary;
                }
            }
        }
        else if(!sortingByX) {   //sortingByY
            if (y<p.get_y()) {
                temporary = *this;
                return temporary;
            }
            else if (y>p.get_y()) {
                temporary = p;
                return temporary;
            }
            else { // point1.y == point2.y
                if (x<p.get_x()) {
                    temporary = *this;
                    return temporary;
                }
                else {
                    temporary = p;
                    return temporary;
                }
            }
        }
        return temporary;
    }
}
