cmake_minimum_required(VERSION 3.5)
project(MEX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(WIN32 AND NOT (${CMAKE_SIZEOF_VOID_P} EQUAL 8))
    message(FATAL_ERROR "You've selected the 32bit version of ${CMAKE_GENERATOR}. \n Please delete the cache (file->Delete Cache) and use the 64bit version. (${CMAKE_GENERATOR} Win64)")
ENDIF()

SET(MATLAB_ROOT $ENV{MATLAB_ROOT} CACHE FILEPATH "directory")

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_definitions(/DMATLAB_MEX_FILE)
add_definitions(/DMX_COMPAT_32)

find_package(Matlab REQUIRED)

IF(MATLAB_FOUND)
    message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
	
	find_package(ZED 3 REQUIRED)
    find_package(CUDA ${ZED_CUDA_VERSION} REQUIRED)
    find_package(OpenCV REQUIRED)
	
    include_directories(include)
	
    add_subdirectory(mex)
ELSE(MATLAB_FOUND)
	message( FATAL_ERROR "MATLAB not found...nothing will be built." )	
ENDIF(MATLAB_FOUND)
