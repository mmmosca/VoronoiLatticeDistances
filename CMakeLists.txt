#[[
Copyright (C) 2018-2022,2026 Marco M. Mosca

This file is part of VoronoiLatticeDistances.

VoronoiLatticeDistances is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VoronoiLatticeDistances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VoronoiLatticeDistances. If not, see <https://www.gnu.org/licenses/>.
*]]
cmake_minimum_required (VERSION 3.8)

project (VoronoiLatticeDistances VERSION 0.8.0 LANGUAGES CXX)

message(STATUS "Project directory: ${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

# set <Package>_DIR env variable to find <Package>Config.cmake files
find_package(Gemmi MODULE REQUIRED)
find_package(Eigen3 MODULE REQUIRED)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/External/boost-1.81.0/install-dir")
file(GLOB BOOST_LIBRARIES
    "${CMAKE_CURRENT_SOURCE_DIR}/External/boost-1.81.0/install-dir/lib/*.lib"
)
find_package(Boost 1.81.0 REQUIRED)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/External/VTK-8.2.0/install-dir")
find_package(VTK 8.2.0 REQUIRED COMPONENTS
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonColor
    vtkCommonTransforms
    vtkFiltersGeneral
    vtkFiltersSources
    vtkIOXML
    vtkInteractionStyle
    vtkRenderingCore
    vtkRenderingOpenGL2
)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/External/CGAL-4.14.3")
find_package(CGAL REQUIRED COMPONENTS Core)

set (WORKDIR "${CMAKE_BINARY_DIR}/Workdir"
		CACHE INTERNAL "${PROJECT_NAME}: WorkDir" FORCE)

add_subdirectory ("CrystaLibrary")
add_subdirectory ("Main")

if (BUILD_TESTING)
    enable_testing()
    include(GoogleTest)
    add_subdirectory ("Tests")
endif()
