cmake_minimum_required(VERSION 3.13)
project(traveling_salesman)

set(CMAKE_CXX_STANDARD 14)

add_executable(traveling_salesman main.cpp EDGE.h POINT.h variant_A.cpp variant_B.cpp)