cmake_minimum_required(VERSION 3.19.0 FATAL_ERROR)
project(cBFGS VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Release)

option(BUILD_SHARED_LIBS "build shared libraries" OFF)

if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not supported. Please create a new directory (e.g. 'Builds') and run CMake from there")
endif ()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(cbfgs)
add_subdirectory(tests)

add_custom_target(__prebuild__ ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_BINARY_DIR}/compile_commands.json" "${PROJECT_SOURCE_DIR}/compile_commands.json"
  )

# vim:fenc=utf-8 et sw=2:
