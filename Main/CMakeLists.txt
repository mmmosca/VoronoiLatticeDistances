# CMakeList.txt: progetto CMake per CrystaLibrary. Includere l'origine e definire
# qui la logica specifica del progetto.
#
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