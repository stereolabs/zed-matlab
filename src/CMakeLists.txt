cmake_minimum_required(VERSION 2.8)
project(MEX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_definitions(/DMATLAB_MEX_FILE)
add_definitions(/DMX_COMPAT_32)

find_package(Matlab REQUIRED)

MESSAGE( STATUS "Matlab path $ENV{MATLAB_ROOT}")

IF(WIN32)	
    SET(ZED_INCLUDE_DIRS $ENV{ZED_INCLUDE_DIRS})
        if (CMAKE_CL_64) # 64 bits
                SET(ZED_LIBRARIES $ENV{ZED_LIBRARIES_64})
        else(CMAKE_CL_64) # 32 bits
                message("32bits compilation is no more available with CUDA7.0")
        endif(CMAKE_CL_64)
    SET(ZED_LIBRARY_DIR $ENV{ZED_LIBRARY_DIR})
    SET(OPENCV_DIR $ENV{OPENCV_DIR})
    SET(Eigen_INCLUDE_DIRS $ENV{EIGEN_DIR})
ELSE()
	find_package(ZED 1.0 REQUIRED)
    SET(Eigen_INCLUDE_DIRS "/usr/include/eigen3") ##depends if installed by packet or by own build install
    SET(Eigen_INCLUDE_DIRS_2 "/usr/local/include/eigen3") ##depends if installed by packet or by own build install
    include_directories(${Eigen_INCLUDE_DIRS_2}) 
ENDIF(WIN32)

find_package(CUDA 7.5 REQUIRED)
find_package(OpenCV 3.1 COMPONENTS core highgui imgproc REQUIRED)

include_directories(${Eigen_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(include)

IF(MATLAB_FOUND)
    message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
    add_subdirectory(mex)
ELSE(MATLAB_FOUND)
    MESSAGE("MATLAB not found...nothing will be built.")
ENDIF(MATLAB_FOUND)
