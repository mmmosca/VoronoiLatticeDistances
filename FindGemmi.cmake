# The following variable should be added to the system variables 
IF(NOT DEFINED ENV{GEMMI_ROOT})
	message(STATUS "Using local package for Gemmi library...")
	set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/External/Gemmi")	
ELSE()
	message(STATUS "Using system environmental variable package for Gemmi library...")
	set(PACKAGE_DIR "$ENV{GEMMI_ROOT}")

ENDIF()

set(GEMMI_INCLUDE_DIR "${PACKAGE_DIR}/include")
set(3RD_INCLUDE_DIR "${PACKAGE_DIR}/third_party" )
set(GEMMI_SRC_DIR "${PACKAGE_DIR}/src")

IF( (NOT EXISTS ${GEMMI_INCLUDE_DIR}) OR (NOT EXISTS ${GEMMI_SRC_DIR}) )
	message(SEND_ERROR "GEMMI directories do not exist, Please assure they are inside the package folder")
	return()
ELSE()
	message(STATUS "Gemmi library detected at: ${PACKAGE_DIR}")
ENDIF()

include_directories(${GEMMI_INCLUDE_DIR} ${3RD_INCLUDE_DIR} )

set (Gemmi_FOUND 1)