# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.8)

project ("VoronoiLatticeDistances")

message(STATUS "Project directory: ${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

# set <Package>_DIR env variable to find <Package>Config.cmake files
find_package(Gemmi REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(VTK REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core MODULE)

set (WORKDIR "${CMAKE_BINARY_DIR}/Workdir"
		CACHE INTERNAL "${PROJECT_NAME}: WorkDir" FORCE)

# Includere i sottoprogetti.
add_subdirectory ("CrystaLibrary")
add_subdirectory ("Main")

