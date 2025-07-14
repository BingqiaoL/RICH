cmake_minimum_required(VERSION 3.10)
project(NegativeWeightGraph)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Include directories
include_directories(${PROJECT_SOURCE_DIR})

# Source files

set(DP_WITH_FILTER_FULL_SOURCES
    dp_with_filter_full.cpp
    directed_graph.cpp
)

set(DP_WITH_FILTER_FULL_SOURCES_ATMOSTK
    dp_with_filter_full_atmost_k.cpp
    directed_graph.cpp
)

# Add executable
add_executable(DPWithFilterFull ${DP_WITH_FILTER_FULL_SOURCES})
add_executable(DPWithFilterFullAtMostk ${DP_WITH_FILTER_FULL_SOURCES_ATMOSTK})
