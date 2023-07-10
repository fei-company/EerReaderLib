cmake_minimum_required(VERSION 3.0)
project(EERReaderLib)

# Set the library version
set(LIBRARY_VERSION_MAJOR 0)
set(LIBRARY_VERSION_MINOR 0)
set(LIBRARY_VERSION_PATCH 1)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Add your source files
set(SOURCES
	./Bitmap.cpp
	./DefectMaskingEER.cpp
	./EerFile.cpp
	./EerFrame.cpp
	./ElectronCountedFramesDecompressor.cpp
	./GainDefectCorrect.cpp
	./stdafx.cpp
)

# Add your header files
set(HEADERS
	./Bitmap.h
	./DefectMaskingEER.h
	./EerFile.h
	./EerFrame.h
	./ElectronCountedFramesDecompressor.h
	./FeiBitStreamer.h
	./GainDefectCorrect.h
	./stdafx.h
	./targetver.h
)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add library target
add_library(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Set library version
set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${LIBRARY_VERSION_MAJOR}.${LIBRARY_VERSION_MINOR}.${LIBRARY_VERSION_PATCH}
    SOVERSION ${LIBRARY_VERSION_MAJOR}
)

# Include directories
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
	${LIBTIFF_INCLUDE_DIR} # Include directory for libtiff
)

# Link directories
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib # Default library directory
    ${LIBTIFF_LIB_DIR} # Library directory for libtiff
)