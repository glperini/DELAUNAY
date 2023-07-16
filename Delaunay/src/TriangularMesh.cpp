#include "TriangularMesh.hpp"

using namespace ProjectLibrary;


namespace ProjectLibrary{

    //getters & setters
    unsigned int TriangularMesh::get_numberOfPoints(){
        return numberOfPoints;
    }

    vector<Point> TriangularMesh::get_points(){
        return points;
    }

    const string TriangularMesh::get_fileName(){
        return fileName;
    }

    const string TriangularMesh::get_outputFilePath(){
        return outputFilePath;
    }

    vector<Triangle> TriangularMesh::get_triangles(){
        return triangles;
    }

    vector<Segment> TriangularMesh::get_segments(){
        return segments;
    }

    Segment TriangularMesh::get_segment(int index){
        return segments[index];
    }

    void TriangularMesh::set_numberOfPoints(unsigned int numberOfPoints){
        this -> numberOfPoints = numberOfPoints;
    }

    void TriangularMesh::set_points(vector<Point> points){
        this -> points = points;
    }
    
    void TriangularMesh::set_fileName(const string fileName){
        this -> fileName = fileName;
    }

    void TriangularMesh::set_outputFilePath(string outputFilePath){
        this -> outputFilePath = outputFilePath;
    }

    void TriangularMesh::set_triangles(vector<Triangle> triangles){
        this -> triangles = triangles;
    }

    void TriangularMesh::set_segments(vector<Segment> segments){
        this -> segments = segments;
    }

    void TriangularMesh::set_segment(Segment segment, int index){
        this -> segments[index] = segment;
    }
    

    // metodi
    bool TriangularMesh::importPoints() {
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

     numberOfPoints = listLines.size();

     if (numberOfPoints == 0)
     {
         cerr << "There are no points" << endl;
         return false;
     }

     points.reserve(numberOfPoints);

     for (const string& line : listLines)
     {
         istringstream converter(line);

        Point p;
        double x;
        double y;
        int id;

         converter >> id >> x >> y;
         p = Point(id, x, y);

         points.push_back(p);

     }

     file.close();
     return true;
    }
    
    void TriangularMesh::drawSegments()
    {
        
        Segment s1 = Segment({triangles[0].get_point(0), triangles[0].get_point(1)});
        Segment s2 = Segment({triangles[0].get_point(1), triangles[0].get_point(2)});
        Segment s3 = Segment({triangles[0].get_point(2), triangles[0].get_point(0)});
        segments.push_back(s1);
        segments.push_back(s2);
        segments.push_back(s3);

        int trianglesSize = triangles.size();
        for (int i=1; i<trianglesSize; i++) {
            //mi creo i tre segmenti del triangolo
            s1.set_points({triangles[i].get_point(0), triangles[i].get_point(1)});
            s2.set_points({triangles[i].get_point(1), triangles[i].get_point(2)});
            s3.set_points({triangles[i].get_point(2), triangles[i].get_point(0)});

            //verifico se è già presente in segments
            bool alreadyExisting1 = false;
            int drawSegmentsSize = segments.size();
            for (int k=0; k<drawSegmentsSize; k++) {
                if (s1.areSameSegment(segments[k])) {
                        alreadyExisting1 = true;
                        break;
                }
            }
            //se non è presente, lo aggiungo
            if (alreadyExisting1 == false)
                segments.push_back(s1);

            bool alreadyExisting2 = false;
            drawSegmentsSize = segments.size();
            for (int k=0; k<drawSegmentsSize; k++) {
                if (s2.areSameSegment(segments[k])) {
                        alreadyExisting2 = true;
                        break;
                }
            }
            if (alreadyExisting2 == false)
                segments.push_back(s2);

            bool alreadyExisting3 = false;
            drawSegmentsSize = segments.size();
            for (int k=0; k<drawSegmentsSize; k++) {
                if (s3.areSameSegment(segments[k])) {
                        alreadyExisting3 = true;
                        break;
                }
            }
            if (alreadyExisting3 == false)
                segments.push_back(s3);
        }
    }

    //Esporta i dati da file
    bool TriangularMesh::exportResult()
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
            file << to_string(segments[i].get_point(0).get_x()) << " "
                 << to_string(segments[i].get_point(0).get_y()) << " "
                 << to_string(segments[i].get_point(1).get_x()) << " "
                 << to_string(segments[i].get_point(1).get_y()) << endl;
        }

        /// Close File
        file.close();

        return true;
    }
}
