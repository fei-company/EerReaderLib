cmake_minimum_required(VERSION 3.0)
project(TestEerReader)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# set output executable directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set ( PROJECT_LINK_LIBS EERReaderLib.lib libtiff_static.lib)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib/Debug ${LIBTIFF_LIB_DIR}) 

include_directories(.. ${LIBTIFF_INCLUDE_DIR}) 

# Sources definition
#set(SOURCES TestEerReader.cpp)

# Target definition
add_executable(TestEerReader TestEerReader.cpp)
target_link_libraries(TestEerReader ${PROJECT_LINK_LIBS} )