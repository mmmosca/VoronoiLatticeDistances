#[[
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*]]
cmake_minimum_required (VERSION 3.8)

project (VoronoiLatticeDistances VERSION 0.7.0 LANGUAGES CXX)

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
