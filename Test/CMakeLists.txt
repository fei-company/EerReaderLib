project(TestEerReader)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib ${LIBTIFF_LIB_DIR})
set (PROJECT_LINK_LIBS EERReaderLib.lib libtiff_static.lib)
include_directories(${LIBTIFF_INCLUDE_DIR}) 

# Target definition
add_executable(TestEerReader TestEerReader.cpp)
target_link_libraries(TestEerReader ${PROJECT_LINK_LIBS} )