cmake_minimum_required(VERSION 3.10)
project(tdsEngLib C)

# ==============================================================================
# 1. OUTPUT DIRECTORIES (MUST BE DEFINED BEFORE 'add_library')
# ==============================================================================
# This ensures that when the library is created, it knows where to go.
# Windows (Visual Studio) will still append a /Release or /Debug folder 
# inside these (e.g., build/bin/Release/tdsEngLib.dll).

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ==============================================================================
# 2. ARCHITECTURE SAFETY CHECK
# ==============================================================================
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR 
        "Error: You are trying to build a 64-bit binary.\n"
        "This library requires 32-bit (x86) to maintain binary compatibility.\n"
        "Please delete the 'build' directory and regenerate using:\n"
        "   cmake -A Win32 ..\n"
    )
endif()

# ==============================================================================
# 3. CONFIGURATION
# ==============================================================================
include_directories(include)

set(SOURCES
    src/tdsEngLib.c
)

# ==============================================================================
# 4. WINDOWS / MSVC SPECIFICS
# ==============================================================================
if(WIN32)
    list(APPEND SOURCES src/tdsEngLib.DEF)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DWIN32_LEAN_AND_MEAN)

    if(MSVC)
        add_compile_options("/arch:IA32")
    endif()
endif()

# ==============================================================================
# 5. BUILD TARGET
# ==============================================================================
# Now that output directories are set, we create the library
add_library(tdsEngLib SHARED ${SOURCES})