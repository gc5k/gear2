cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -pthread -fopenmp")
set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_BUILD_TYPE "Debug")

project(proEigenGWAS)

include_directories(include)
#cmake include boost
find_package(Boost 1.72.0 COMPONENTS program_options REQUIRED)
find_package(Eigen3 3 REQUIRED NO_MODULE)

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})

add_executable(gear src/gear.cpp src/genotype.cpp)
target_include_directories(gear PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(gear ${Boost_LIBRARIES} Eigen3::Eigen)
