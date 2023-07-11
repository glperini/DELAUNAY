#include "empty_class.hpp"
//#include "empty_class.cpp"

using namespace ProjectLibrary;

//bool importPoints(TriangularMesh& mesh, const string& fileName);

int main()
{
    TriangularMesh mesh;

    //importo i punti
    if(!mesh.importPoints(mesh, "./Test2.csv"))
        return 1;

    vector<Point> points = mesh.points;

    //ordino i punti per x crescente e y crescente secondo l'algoritmo Mergesort O(nlogn)
    MergeSorter toBeSortedX;
    MergeSorter toBeSortedY;
    int sx = 0;
    int dx = points.size()-1;
    bool sortingByX = true;
    vector<Point> sortedX = toBeSortedX.MergeSort(points, sx, dx, sortingByX);
    sortingByX = !sortingByX;
    vector<Point> sortedY = toBeSortedY.MergeSort(points, sx, dx, sortingByX);

    //inizio l'algoritmo creando il primo triangolo
    vector<Triangle> triangles;
    vector<Segment> convexHull;
    vector<int> hullTrianglesIndices;
    vector<Point> hullPoints;
    Delaunay delaunay;

    //per creare il primo triangolo
    delaunay.firstTriangle(sortedX,
                           sortedY,
                           triangles,
                           convexHull,
                           hullPoints,
                           hullTrianglesIndices);


    //inserisco un nuovo punto nella triangolazione
    for (unsigned int p=0; p<sortedX.size(); p++) {
        int triangleIndex = 0;
        Point newPoint = sortedX[p];
        if (delaunay.findTriangleContainingPoint(newPoint, triangles, triangleIndex) != -1) {

            //il nuovo punto inserito si trova all'interno della triangolazione
            delaunay.inTriangle(newPoint,
                                triangles,
                                triangleIndex,
                                hullTrianglesIndices);
        }
        else {

            //il nuovo punto inserito si trova all'esterno della triangolazione
            delaunay.outTriangle(newPoint,
                                 triangles,
                                 convexHull,
                                 hullPoints,
                                 hullTrianglesIndices);
        }
    }

    vector<Segment> segments = delaunay.drawSegments(triangles);

    if (!mesh.exportResult("./Results.csv", segments))
        return 1;

    return 0;

}
