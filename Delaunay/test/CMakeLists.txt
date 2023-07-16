list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Point.hpp)
list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Segment.hpp)
list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Triangle.hpp)
list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Delaunay.hpp)
list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_TriangularMesh.hpp)
list(APPEND delaunay_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_HeapSort.hpp)

list(APPEND delaunay_test_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(delaunay_test_headers ${delaunay_test_headers} PARENT_SCOPE)
set(delaunay_test_includes ${delaunay_test_includes} PARENT_SCOPE)

