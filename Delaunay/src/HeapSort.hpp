#ifndef __HEAPSORT_H
#define __HEAPSORT_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include <cmath>
#include <algorithm>

#include "Point.hpp"

using namespace std;
using namespace Eigen;

using namespace ProjectLibrary;

namespace SortingLibrary {

void heapify(vector<Point>& v, int n, int i, bool sortingByX)
{
    int parentIndex = i;
    int leftChildIndex = 2*i + 1;
    int rightChildIndex = 2*i + 2;

    if (leftChildIndex < n && v[leftChildIndex].getOrder(v[parentIndex], sortingByX).get_x() == v[parentIndex].get_x()
            && v[leftChildIndex].getOrder(v[parentIndex], sortingByX).get_y() == v[parentIndex].get_y())
        parentIndex = leftChildIndex;

    if (rightChildIndex < n && v[rightChildIndex].getOrder(v[parentIndex], sortingByX).get_x() == v[parentIndex].get_x()
            && v[rightChildIndex].getOrder(v[parentIndex], sortingByX).get_y() == v[parentIndex].get_y() )
        parentIndex = rightChildIndex;

    if (parentIndex!=i) {
        swap(v[i], v[parentIndex]);
        heapify(v, n, parentIndex, sortingByX);
    }
}

vector<Point> heapSort(vector<Point> v, bool sortingByX)
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
}

#endif // __HEAPSORT_H
