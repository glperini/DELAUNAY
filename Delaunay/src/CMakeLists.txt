list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Point.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Segment.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Triangle.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/TriangularMesh.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Delaunay.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/HeapSort.hpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Point.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Segment.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Triangle.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Triangle.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/TriangularMesh.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/Delaunay.cpp)

list(APPEND delaunay_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(delaunay_sources ${delaunay_sources} PARENT_SCOPE)
set(delaunay_headers ${delaunay_headers} PARENT_SCOPE)
set(delaunay_includes ${delaunay_includes} PARENT_SCOPE)
