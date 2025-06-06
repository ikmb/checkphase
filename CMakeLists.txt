cmake_minimum_required(VERSION 3.9)

# Import local CMakeModules folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

project(checkphase C CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED 17)

#############################################################################
# C++ source files and target definition
#############################################################################
set(cpp_sources
    checkphase.cpp
    version.cpp
    Args.cpp
)

#############################################################################
# Set up Git version extraction
#############################################################################
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.template
    ${CMAKE_BINARY_DIR}/generated/version.h
    )

include_directories("${CMAKE_BINARY_DIR}/generated")

#############################################################################
# Configure libraries
#############################################################################

add_definitions(-x c++)
add_executable(${PROJECT_NAME} ${cpp_sources})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

# Zlib
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

# HTSlib
find_package(HTSlib REQUIRED)
include_directories(${HTSlib_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSlib_LIBRARIES})

# Boost
find_package(Boost 1.65.1 REQUIRED
    COMPONENTS
    filesystem
    program_options
)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

target_link_libraries(${PROJECT_NAME} curl)

# additional CXX options
set(CXX_WARNINGS "-Wall -Wextra -Wshadow -Wno-deprecated-declarations -Wmissing-include-dirs -Wstrict-aliasing=2 -Wfloat-equal -Wcast-align")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_OPTIONS} ${cxx_compiler_opts} -std=c++11 ${CXX_WARNINGS} -fopenmp -march=native")
