#include "Delaunay.hpp"
#include "TriangularMesh.hpp"
#include "HeapSort.hpp"

using namespace ProjectLibrary;
using namespace SortingLibrary;

int main()
{
    TriangularMesh mesh;

    //importo i punti
    mesh.set_fileName("./Test2.csv");
    if(!mesh.importPoints())
        return 1;

    //ordino i punti per x crescente e y crescente secondo l'algoritmo Heapsort O(nlogn)
    bool sortingByX = true;
    vector<Point> sortedX = heapSort(mesh.get_points(), sortingByX);
    sortingByX = !sortingByX;
    vector<Point> sortedY = heapSort(mesh.get_points(), sortingByX);

    //inizio l'algoritmo creando il primo triangolo
    Delaunay delaunay;

    //per creare il primo triangolo
    delaunay.firstTriangle(sortedX,
                           sortedY);

    mesh.set_triangles(delaunay.get_triangles());

    //inserisco un nuovo punto nella triangolazione
    for (unsigned int p=0; p<sortedX.size(); p++) {
        delaunay.set_triangleIndex(0);
        Point newPoint = sortedX[p];
        if (delaunay.findTriangleContainingPoint(newPoint) != -1) {

            //il nuovo punto inserito si trova all'interno della triangolazione
            delaunay.inTriangle(newPoint);
        }
        else {

            //il nuovo punto inserito si trova all'esterno della triangolazione
            delaunay.outTriangle(newPoint);
        }

        mesh.set_triangles(delaunay.get_triangles());

    }

    mesh.drawSegments();
    mesh.set_outputFilePath("Results.csv");
    if (!mesh.exportResult())
        return 1;

    return 0;
}
