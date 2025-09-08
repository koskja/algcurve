cmake_minimum_required(VERSION 3.10)

project(algcurve LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(INCLUDE_DEBUG_SYMBOLS "Include debug symbols")
option(ENABLE_SANITIZERS "Enable Address and Undefined Behavior Sanitizers")
option(FRAME_POINTERS "Enable frame pointers")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(ENABLE_SANITIZERS ON)
    set(INCLUDE_DEBUG_SYMBOLS ON)
endif()
if(ENABLE_SANITIZERS OR INCLUDE_DEBUG_SYMBOLS)
    set(FRAME_POINTERS ON)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC
  main.cpp
  input.cpp
  render.cpp
  thread.cpp
  paper.cpp
)

add_executable(algcurve ${SRC})
target_include_directories(algcurve PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(algcurve PRIVATE
    -Wall -Wextra -Wshadow -Wconversion
  )
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(algcurve PRIVATE -O3)
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(algcurve PRIVATE -O0)
  endif()
  if(INCLUDE_DEBUG_SYMBOLS)
    target_compile_options(algcurve PRIVATE -g3)
  endif()
  if(FRAME_POINTERS)
      target_compile_options(algcurve PRIVATE -fno-omit-frame-pointer)
  endif()
  if(ENABLE_SANITIZERS)
    target_compile_options(algcurve PRIVATE -fsanitize=address,undefined)
    target_link_options(algcurve PRIVATE -fsanitize=address,undefined)
  endif()
  
  if(NOT WIN32)
    target_compile_options(algcurve PRIVATE -march=native)
  endif()
elseif(MSVC)
  target_compile_options(algcurve PRIVATE /W4)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(algcurve PRIVATE /Ox)
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(algcurve PRIVATE /Od)
  endif()
  if(INCLUDE_DEBUG_SYMBOLS)
    target_compile_options(algcurve PRIVATE /Zi)
  endif()
  if(ENABLE_SANITIZERS)
    target_compile_options(algcurve PRIVATE /fsanitize=address)
    target_link_options(algcurve PRIVATE /fsanitize=address)
  endif()
endif()

if(NOT WIN32)
  target_link_libraries(algcurve PRIVATE m)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


