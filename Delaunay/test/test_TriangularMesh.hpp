#ifndef __TEST_TRIANGULARMESH_H
#define __TEST_TRIANGULARMESH_H

#include <gtest/gtest.h>
#include "TriangularMesh.hpp"

using namespace ProjectLibrary;

using namespace testing;

///TEST TriangularMesh
TEST(TestTriangularMesh, Test_drawSegments) {

    TriangularMesh mesh;

    Point p1 = Point({3.0, 3.0});
    Point p2 = Point({5.0, 7.0});
    Point p3 = Point({11.0, 3.0});
    Point p4 = Point({6.0, 2.0});
    Point p5 = Point({3.0, 6.0});

    Triangle triangle1 = Triangle({p1,p3,p2});
    Triangle triangle2 = Triangle({p1,p4,p3});
    Triangle triangle3 = Triangle({p2,p5,p1});

    mesh.set_triangles({triangle1, triangle2, triangle3});

    mesh.drawSegments();

    vector<Segment> expectedSegments = {Segment({p1,p3}), Segment({p3,p2}), Segment({p2,p1}), Segment({p1,p4}), Segment({p4,p3}), Segment({p2,p5}), Segment({p5,p1})};

    for (int i = 0; i<7; i++)
        EXPECT_EQ(mesh.get_segment(i).areSameSegment(expectedSegments[i]), true);
}


#endif // __TEST_TRIANGULARMESH_H
