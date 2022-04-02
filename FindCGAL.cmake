# The following variable should be added to the system variables 

IF(NOT DEFINED ENV{CGAL_DIR})
	message(STATUS "Using local package for CGAL library...")
	set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/External/CGAL-4.14.3")
	find_package(CGAL REQUIRED PATHS "${PACKAGE_DIR}" COMPONENTS Core CONFIG)
	message(STATUS "CGAL library detected at: ${PACKAGE_DIR}")
ELSE()
	message(STATUS "Using system environmental variable package for CGAL library...")
	find_package(CGAL REQUIRED COMPONENTS Core CONFIG)
ENDIF()
