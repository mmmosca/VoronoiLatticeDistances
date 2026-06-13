#[[Copyright (C) 2018-2022,2026 Marco M. Mosca

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
]]
project(CrystaLibrary)
cmake_minimum_required (VERSION 3.8)

set( CURR_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
set( INCLUDE_DIR "${CURR_DIR}/include" )
set( SRC_DIR "${CURR_DIR}/src" )

# Aggiungere l'origine all'eseguibile di questo progetto.

IF(WIN32)
	# Add 'lib' at the beginning of the file
	# set(CMAKE_SHARED_LIBRARY_PREFIX "lib")

	# Export symbols in a .lib library file
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
ENDIF()

include(${CGAL_USE_FILE})
include(${VTK_USE_FILE})

file(GLOB CFILES ${SRC_DIR}/*.cpp)
file(GLOB HFILES ${INCLUDE_DIR}/*.h)

if (SHARED)
	add_library(${PROJECT_NAME} SHARED ${CFILES} ${HFILES})
else()
	add_library(${PROJECT_NAME} STATIC ${CFILES} ${HFILES})
ENDIF()

target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR} )
target_link_libraries(${PROJECT_NAME} ${BOOST_LIBRARIES} ${CGAL_LIBRARIES} ${VTK_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${WORKDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${WORKDIR}"
	ARCHIVE_OUTPUT_DIRECTORY "${WORKDIR}"
)

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
install(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION lib)
install(FILES ${HFILES} DESTINATION include)
