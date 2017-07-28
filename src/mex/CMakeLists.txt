set(CPP_FILE mexZED)

include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${ZED_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

link_directories(${ZED_LIBRARY_DIR})
link_directories(${OpenCV_LIBRARY_DIRS})
link_directories(${CUDA_LIBRARY_DIRS})

set(CPP_FILE_NAME ${CPP_FILE}.cpp)
INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
add_library(${CPP_FILE} SHARED ${CPP_FILE_NAME} ${CMAKE_SOURCE_DIR}/Matlabdef.def)

add_definitions(-std=c++0x)

target_link_libraries(${CPP_FILE}
					${MATLAB_LIBRARIES}
					${ZED_LIBRARIES}
					${OpenCV_LIBRARIES}
					${CUDA_LIBRARIES} ${CUDA_nppi_LIBRARY} ${CUDA_npps_LIBRARY})

if(WIN32)
    SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexw64)
else(WIN32)
    if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexa64 PREFIX "")
    else(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexglx PREFIX "")
    endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
endif(WIN32)

install(TARGETS ${CPP_FILE} DESTINATION ../matlab)
