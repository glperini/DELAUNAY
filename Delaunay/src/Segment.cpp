#include "Segment.hpp"

namespace ProjectLibrary {

    //Costruttore p1,p2
    Segment::Segment(vector<Point> points):
        points(points)
    {}

    //getters & setters
    vector<Point> Segment::get_points(){
        return points;
    }

    Point Segment::get_point(int index){
        return points[index];
    }

    Vector2d Segment::get_coefficients(){
        return coefficients;
    }

    double Segment::get_coefficient(int index){
        return coefficients[index];
    }

    void Segment::set_points(vector<Point> points){  
        this -> points = points;
    }

    void Segment::set_point(Point point, int index){
        this -> points[index] = point;
    }   

    void Segment::set_coefficients(Vector2d coefficients){
        this -> coefficients = coefficients;
    }

    void Segment::set_coefficient(double coefficient, int index){
        this -> coefficients[index] = coefficient;
    }

    //metodi 
    void Segment::lineCoefficients(){
        //se i due punti non sono allineati verticalmente (non hanno la stessa x):
        if (get_point(0).get_x()!=get_point(1).get_x()){
            //Calcoliamo m
            double m = (get_point(1).get_y() - get_point(0).get_y()) / (get_point(1).get_x() - get_point(0).get_x());
            //Calcoliamo q
            double q = get_point(0).get_y() - m*get_point(0).get_x();

            //calcolati m e q, mettiamoli dentro il vettore coefficients
            coefficients[0] = m;
            coefficients[1] = q;
        }
        else {  //se hanno la stessa x
            coefficients[0] = numeric_limits<double>::quiet_NaN();
            coefficients[1] = get_point(0).get_x();  //in questo caso coefficients[1] è l'intersezione con l'asse delle x
        }
    }

    double Segment::getAngle(Segment segment){
        double angleDegree;
        lineCoefficients();
        segment.lineCoefficients();
        Vector2d coefficients2 = segment.get_coefficients();
         //verifico se i due segmenti hanno stesso m (cioè sono paralleli) oppure se sono entrambi verticali
        if (coefficients[0]==coefficients2[0] || (isnan(coefficients[0]) && isnan(coefficients2[0])))
            angleDegree = 0.0;
        //se i due segmenti sono perpendicolari (ossia m1*m2 = -1):
        else if (coefficients[0]*coefficients2[0]==-1)
            angleDegree = 90.0;
        //se i due segmenti non sono nè paralleli nè perpendicolari:
        else {
            double vector1X = get_point(1).get_x() - get_point(0).get_x();
            double vector1Y = get_point(1).get_y() - get_point(0).get_y();
            double vector2X = segment.get_point(1).get_x() - segment.get_point(0).get_x();
            double vector2Y = segment.get_point(1).get_y() - segment.get_point(0).get_y();

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

    bool Segment::areSameSegment(Segment segment){
        //vector<Point> points2 = segment.get_points()
        if((points[0].areSamePoint(segment.get_point(0)) && points[1].areSamePoint(segment.get_point(1))) ||
           (points[0].areSamePoint(segment.get_point(1)) && points[1].areSamePoint(segment.get_point(0)))) {
            return true;
        }
        else
            return false;
    }
}
