cmake_minimum_required(VERSION 3.10)
project(tdsEngLib C)

# 1. Setup Include Directories
# Allows code in 'src' to find headers in 'include'
include_directories(include)

# 2. Define Source Files
set(SOURCES
    src/tdsEngLib.c
)

# 3. Handle Windows DLL Exports (.DEF file)
if(WIN32)
    # Add the definition file to sources for Windows builds
    list(APPEND SOURCES src/tdsEngLib.DEF)
    
    # Standard Windows definitions
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DWIN32_LEAN_AND_MEAN)
endif()

# 4. Create the Shared Library (DLL)
add_library(tdsEngLib SHARED ${SOURCES})

# 5. Output Directory Settings (Optional but helpful)
# Puts the .dll and .lib in a 'bin' and 'lib' folder respectively
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)