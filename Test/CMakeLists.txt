project(TestEerReader)

if(WIN32)
    set (PROJECT_LINK_LIBS EERReaderLib.lib libtiff_static.lib)
else()
    set (PROJECT_LINK_LIBS libEERReaderLib.a libtiff.a)
endif()

include_directories(${LIBTIFF_INCLUDE_DIR})

# Target definition
add_executable(TestEerReader TestEerReader.cpp)
target_link_libraries(TestEerReader ${PROJECT_LINK_LIBS} )