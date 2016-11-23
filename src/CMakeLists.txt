cmake_minimum_required(VERSION 2.8)
project(MEX)

IF(WIN32 AND NOT CMAKE_GENERATOR MATCHES "[a-z A-Z 0-9]*64" )
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
	
    IF(WIN32)
        SET(ZED_INCLUDE_DIRS $ENV{ZED_INCLUDE_DIRS})
        SET(ZED_LIBRARIES $ENV{ZED_LIBRARIES_64})
        SET(ZED_LIBRARY_DIR $ENV{ZED_LIBRARY_DIR})
        SET(OPENCV_DIR $ENV{OPENCV_DIR})
        SET(EIGEN3_INCLUDE_DIR $ENV{EIGEN_DIR})
    ELSE()
        find_package(ZED 1.0 REQUIRED)
        find_package(Eigen3 REQUIRED)
    ENDIF(WIN32)

    find_package(CUDA 8 REQUIRED)
    find_package(OpenCV 3.1 COMPONENTS core highgui imgproc REQUIRED)

    include_directories(${EIGEN3_INCLUDE_DIR})
    include_directories(include)

    add_subdirectory(mex)
ELSE(MATLAB_FOUND)
	message( FATAL_ERROR "MATLAB not found...nothing will be built." )	
ENDIF(MATLAB_FOUND)
