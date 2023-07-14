#include "objects.hpp"

namespace LibraryObjects
{

///DEF Triangle
//Costruttore
Triangle::Triangle(Point& p1, Point& p2, Point& p3):
        p1(p1),
        p2(p2),
        p3(p3)
{}

//Calcola area triangolo con formula di Gauss
double Triangle::Area(Point& p1, Point& p2, Point& p3)
{
    double Area = 0;
    Area = p1.x * p2.y - p2.x * p1.y +
           p2.x * p3.y - p3.x * p2.y +
           p3.x * p1.y - p1.x * p3.y;
    Area = 0.5 * abs(Area);
    return Area;
}


///DEF Point
//Costruttore
Point::Point(int& id, double& x, double& y):
  id(id),
  x(x),
  y(y)
{}

//Costruttore
Point::Point(double x, double y):
    x(x),
    y(y)
{id=0;}


///DEF Segment
//Costruttore
Segment::Segment(Point& p1, Point& p2,  Vector2d& coefficients):
   p1(p1),
   p2(p2),
   coefficients(coefficients)
{}

//Costruttore
Segment::Segment(Point p1, Point p2):
   p1(p1),
   p2(p2)
{}

//Calcola i coefficienti di una retta (m e q)
Vector2d Segment::lineCoefficients(Point& p1, Point& p2){
   Vector2d coefficients;

   //se i due punti non sono allineati verticalmente (non hanno la stessa x):
   if (p1.x!=p2.x){
   //Calcoliamo m
   double m = (p2.y - p1.y) / (p2.x - p1.x);
   //Calcoliamo q
   double q = p1.y - m*p1.x;

   //calcolati m e q, mettiamoli dentro il vettore coefficients
   coefficients[0] = m;
   coefficients[1] = q;;
   }
   else {  //se hanno la stessa x
       coefficients[0] = numeric_limits<double>::quiet_NaN();
       coefficients[1] = p1.x;  //in questo caso coefficients[1] è l'intersezione con l'asse delle x
   }

   return coefficients;
 }


 ///DEF TriangularMesh
 //Importa i dati da file
 bool TriangularMesh::importPoints(TriangularMesh& mesh, const string& fileName) {
     //Importo i dati: Id, x, y
     ifstream file;
     file.open(fileName);

     if(file.fail())
         return false;

     list<string> listLines;
     string line;
     while (getline(file, line))
         listLines.push_back(line);

     file.close();

     listLines.pop_front();

     mesh.NumberOfPoints = listLines.size();

     if (mesh.NumberOfPoints == 0)
     {
         cerr << "There are no points" << endl;
         return false;
     }

     mesh.points.reserve(mesh.NumberOfPoints);

     for (const string& line : listLines)
     {
         istringstream converter(line);

         Point coord;

         converter >> coord.id >> coord.x >> coord.y;

         mesh.points.push_back(coord);

     }

     file.close();
     return true;
 }

 //Esporta i dati da file
 bool TriangularMesh::exportResult(const string& outputFilePath, vector<Segment>& segments)
 {
   /// Open File
   ofstream file;
   file.open(outputFilePath);

   if (file.fail())
   {
     cerr<< "file open failed"<< endl;
     return false;
   }

   int segmentsSize = segments.size();
   file << "p1.x p1.y p2.x p2.y" << endl;
   for (int i=0; i<segmentsSize; i++) {
    file << to_string(segments[i].p1.x) << " " << to_string(segments[i].p1.y) << " " <<  to_string(segments[i].p2.x) << " " <<  to_string(segments[i].p2.y) << endl;
   }

   /// Close File
   file.close();

   return true;
 }

}
