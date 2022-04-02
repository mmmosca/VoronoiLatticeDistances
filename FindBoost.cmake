# The following variable should be added to the system variables 
# BOOST_ROOT: Root of BOOST installation

IF(NOT DEFINED ENV{BOOST_ROOT})
	message(SEND_ERROR "BOOST_ROOT not defined, Please add BOOST_ROOT as environmental system variable.")
	return()
ENDIF()

set(BOOST_INCLUDE_DIR "$ENV{BOOST_ROOT}")
set(BOOST_LIB_DIR "$ENV{BOOST_ROOT}/lib32-msvc-14.1")
set(BOOST_BINARY_DIR "$ENV{BOOST_ROOT}/lib32-msvc-14.1")

IF( (NOT EXISTS ${BOOST_INCLUDE_DIR}) OR (NOT EXISTS ${BOOST_LIB_DIR}) )
	message(SEND_ERROR "BOOST directories do not exist, Please assure they are inside the package folder")
	return()
ENDIF()

#FILE(GLOB BOOST_DYN_LIBS ${BOOST_LIB_DIR}/*.lib)
#FILE(GLOB BOOST_DYN_DLLS ${BOOST_BINARY_DIR}/*.dll)


include_directories(${BOOST_INCLUDE_DIR})
link_directories(${BOOST_LIB_DIR})

set (Boost_FOUND 1)