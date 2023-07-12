#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

#include "objects.hpp"

using namespace LibraryObjects;

using namespace std;
using namespace Eigen;

namespace LibraryFunctions
{

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

#endif // __FUNCTIONS_H
