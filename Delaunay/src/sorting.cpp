#include "objects.hpp"
#include "sorting.hpp"

using namespace LibraryObjects;

namespace LibrarySorting
{

///DEF HeapSorter
void HeapSorter::heapify(vector<Point>& v, int n, int i, bool sortingByX)
{
    int parentIndex = i;
    int leftChildIndex = 2*i + 1;
    int rightChildIndex = 2*i + 2;

    if (leftChildIndex < n && getOrder(v[leftChildIndex], v[parentIndex], sortingByX).x == v[parentIndex].x
                           && getOrder(v[leftChildIndex], v[parentIndex], sortingByX).y == v[parentIndex].y)
        parentIndex = leftChildIndex;

    if (rightChildIndex < n && getOrder(v[rightChildIndex], v[parentIndex], sortingByX).x == v[parentIndex].x
                            && getOrder(v[rightChildIndex], v[parentIndex], sortingByX).y == v[parentIndex].y )
        parentIndex = rightChildIndex;

    if (parentIndex!=i) {
        swap(v[i], v[parentIndex]);
        heapify(v, n, parentIndex, sortingByX);
    }
}

vector<Point> HeapSorter::heapSort(vector<Point>& v, bool sortingByX)
{
    int n = v.size();

    //Enqueue
    for (int i=n/2-1; i>=0; i--)        //a partire dall'ultimo nodo padre del heap eseguo heapify fino in cima
        heapify(v, n, i, sortingByX);

    //Dequeue
    for (int i=n-1; i>0; i--) {
        swap(v[0], v[i]);
        heapify(v, i, 0, sortingByX);
    }

    return v;
}

//Definisce ordinamento tra punti
Point HeapSorter::getOrder(Point& point1, Point& point2, bool& sortingByX){
    if (sortingByX == true) {
        if (point1.x<point2.x)
            return point1;
        else if (point1.x>point2.x)
            return point2;
        else { // point1.x == point2.x
            if (point1.y<point2.y)
                return point1;
            else
                return point2;
        }
    }
    else {   //sortingByY
        if (point1.y<point2.y)
            return point1;
        else if (point1.y>point2.y)
            return point2;
        else { // point1.y == point2.y
            if (point1.x<point2.x)
                return point1;
            else
                return point2;
        }
    }
}

}
