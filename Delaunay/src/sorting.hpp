#ifndef __SORTING_H
#define __SORTING_H

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

namespace LibrarySorting
{

//DICHIARAZIONE STRUTTURA HeapSorter
struct HeapSorter
    {
    int n = 0;
    int i = 0;
    bool sortingByX = false;
    vector<Point> v;

    void heapify(vector<Point>& v, int n, int i, bool sortingByX);
    vector<Point> heapSort(vector<Point>& v, bool sortingByX);
    Point getOrder(Point& point1, Point& point2, bool& sortingByX);
    };

}

#endif // __SORTING_H
