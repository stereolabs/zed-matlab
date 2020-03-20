set(TARGET mexZED)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()



include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${ZED_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

link_directories(${ZED_LIBRARY_DIR})
link_directories(${OpenCV_LIBRARY_DIRS})
link_directories(${CUDA_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
add_library(${TARGET} SHARED mexZED.cpp ${CMAKE_SOURCE_DIR}/Matlabdef.def)

add_definitions(-std=c++0x)

target_link_libraries(${TARGET}
					${MATLAB_LIBRARIES}
					${ZED_LIBRARIES}
					${OpenCV_LIBRARIES}
					${CUDA_LIBRARIES} ${CUDA_NPP_LIBRARIES_ZED})

if(WIN32)
    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexw64)
else(WIN32)
    if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexa64 PREFIX "")
    else(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexglx PREFIX "")
    endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
endif(WIN32)

install(TARGETS ${TARGET} DESTINATION ../matlab)
