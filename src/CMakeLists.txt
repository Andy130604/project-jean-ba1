cmake_minimum_required(VERSION 3.5)

project(PATHWAYSEARCH-2023)

# set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-Wall")

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB PROJECT_SOURCES
     "*.cpp"
     "*.hpp"
      )
# Find any version 2.X of SFML, first trying 2.5 or above (for which CMake configuration changed)

add_executable (pathsearch  ${PROJECT_SOURCES})
target_link_libraries(pathsearch)

