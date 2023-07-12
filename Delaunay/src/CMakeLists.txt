list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/objects.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/sorting.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/delaunay.hpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/functions.hpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/objects.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/sorting.cpp)
list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/delaunay.cpp)

list(APPEND delaunay_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(delaunay_sources ${delaunay_sources} PARENT_SCOPE)
set(delaunay_headers ${delaunay_headers} PARENT_SCOPE)
set(delaunay_includes ${delaunay_includes} PARENT_SCOPE)
