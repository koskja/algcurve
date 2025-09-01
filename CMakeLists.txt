cmake_minimum_required(VERSION 3.20)

# Project
project(algcurve LANGUAGES CXX)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Option to include debug symbols in release builds
option(INCLUDE_DEBUG_SYMBOLS "Include debug symbols in release builds" OFF)
set(INCLUDE_DEBUG_SYMBOLS ON)

# Require C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler selection: favor Clang
if(WIN32)
  # On Windows we prefer clang-cl with MSVC-like command line (works with VS IDE)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT MSVC)
    message(STATUS "If you want clang-cl, configure with -T clangcl or use the clang-cl preset.")
  endif()
else()
  # On non-Windows we prefer clang++
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Consider configuring with CXX=clang++ for consistent behavior.")
  endif()
endif()

# Sources
set(SRC
  main.cpp
  input.cpp
  render.cpp
  thread.cpp
  paper.cpp
)

add_executable(algcurve ${SRC})
target_include_directories(algcurve PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Common warnings and opts
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(algcurve PRIVATE
    -Wall -Wextra -Wshadow -Wconversion
    -ffast-math
  )

  # Optimization flags based on build type
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(algcurve PRIVATE -O3)
    if(INCLUDE_DEBUG_SYMBOLS)
      target_compile_options(algcurve PRIVATE -g)
    endif()
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(algcurve PRIVATE -O0 -g)
  endif()
  
  if(NOT WIN32)
    target_compile_options(algcurve PRIVATE -march=native)
  endif()
elseif(MSVC)
  # If building with MSVC (or clang-cl under MSVC), use corresponding flags
  target_compile_options(algcurve PRIVATE /W4)
  
  # Optimization flags based on build type
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(algcurve PRIVATE /O2)
    if(INCLUDE_DEBUG_SYMBOLS)
      target_compile_options(algcurve PRIVATE /Zi)
    endif()
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(algcurve PRIVATE /Od /Zi)
  endif()
endif()

# Defines/compat per platform
if(WIN32)
  # Nothing special for Windows yet; link libraries if needed later
else()
  # Link libm explicitly on Unix-like
  target_link_libraries(algcurve PRIVATE m)
endif()

# Enable compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


