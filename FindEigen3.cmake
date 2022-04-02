# The following variable should be added to the system variables 

IF(NOT DEFINED ENV{EIGEN3_ROOT})
	message(STATUS "Using local package for Eigen3 library...")
	set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/External/eigen-3.3.7")	
ELSE()
	message(STATUS "Using system environmental variable package for Eigen3 library...")
	set(PACKAGE_DIR "$ENV{EIGEN3_ROOT}")

ENDIF()

set(EIGEN3_INCLUDE_DIR "${PACKAGE_DIR}")
set(EIGEN3_SRC_DIR "${PACKAGE_DIR}/Eigen/src")

IF( (NOT EXISTS ${EIGEN3_INCLUDE_DIR}) OR (NOT EXISTS ${EIGEN3_SRC_DIR}) )
	message(SEND_ERROR "EIGEN3 directories do not exist, Please assure they are inside the package folder")
	return()
ELSE()
	message(STATUS "Eigen library detected at: ${PACKAGE_DIR}")
ENDIF()

include_directories(${EIGEN3_INCLUDE_DIR})

set (Eigen3_FOUND 1)