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
project(main)

add_executable(voronoilatticedistances compute_voronoidistances.cpp)
add_executable(voronoilatticedistances_off compute_voronoidistances_fromOFF.cpp)

add_dependencies(voronoilatticedistances CrystaLibrary)
add_dependencies(voronoilatticedistances_off CrystaLibrary)

target_link_libraries(voronoilatticedistances CrystaLibrary)
target_link_libraries(voronoilatticedistances_off CrystaLibrary)


set_target_properties(voronoilatticedistances PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${WORKDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${WORKDIR}"
	ARCHIVE_OUTPUT_DIRECTORY "${WORKDIR}")

set_target_properties(voronoilatticedistances_off PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${WORKDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${WORKDIR}"
	ARCHIVE_OUTPUT_DIRECTORY "${WORKDIR}")

install(TARGETS voronoilatticedistances voronoilatticedistances_off DESTINATION bin)
