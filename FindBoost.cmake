#[[*
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