# GBTL Algorithms in SuiteSparse:GraphBLAS,
# Scott McMillan, (c) 2017, All Rights Reserved.
#
# CMakeLists.tx: instructions for cmake to build GBTL algorithms.
# An ANSI C11 compiler is required.
#
# To compile:
#
#   cd build
#   cmake ..
#   make
#
# If that fails for any reason, make sure your compiler support ANSI C11.  You
# could try changing your compiler, for example, use one of these commands:
#
#   CC=icc cmake ..
#   CC=xlc cmake ..
#   CC=gcc cmake ..
#
# To run the demos:
#
#   cd ../Demo
#   ./<executable file name>
#

# cmake 3.0 is prefered.
cmake_minimum_required ( VERSION 2.8.12 )

project ( HPEC18_Tutorial )

include ( GNUInstallDirs )

if ( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy ( SET CMP0042 NEW )
endif ( )

if (NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release )
endif ( )

# For development only, not for end-users:
# set ( CMAKE_BUILD_TYPE Debug )

set ( CMAKE_INCLUDE_CURRENT_DIR ON )

# include directories for both graphblas and graphblasdemo libraries
include_directories ( include )
link_directories ( lib )
find_library(graphblas libgraphblas)

# check which compiler is being used.  If you need to make
# compiler-specific modifications, here is the place to do it.
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # cmake 2.8 workaround: gcc needs to be told to do ANSI C11.
    # cmake 3.0 doesn't have this problem.
    set (CMAKE_C_FLAGS  "-std=gnu11 -lm -fopenmp")
    if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.9)
        message (FATAL_ERROR "gcc version must be at least 4.9")
    endif ( )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    # options for icc: also needs -std=c11
    set (CMAKE_C_FLAGS  "-std=c11 -fopenmp")
    if (CMAKE_C_COMPILER_VERSION VERSION_LESS 18.0)
        message (FATAL_ERROR "icc version must be at least 18.0")
    endif ( )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    # options for clang
    if (CMAKE_C_COMPILER_VERSION VERSION_LESS 3.3)
        message (FATAL_ERROR "clang version must be at least 3.3")
    endif ( )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    # options for MicroSoft Visual Studio
endif ( )

# create the dynamic GBTL algorithms library.  Requires ANSI C11
#file ( GLOB GBTL_ALGS_SOURCES "GBTL/*.c" "Source/Generated/*.c" )
#add_library ( gbtl_algorithms SHARED ${GBTL_ALGS_SOURCES} )
#SET_TARGET_PROPERTIES ( graphblas PROPERTIES VERSION 1.2.0
#    SOVERSION 1
#    C_STANDARD_REQUIRED 11
#    PUBLIC_HEADER "GBTL/algorithms.h" )
#set_property ( TARGET graphblas PROPERTY C_STANDARD 11 )

# create the static GBTL algorithms library.  Requires ANSI C11
#add_library ( gbtl_algorithms_static STATIC ${GBTL_ALGS_SOURCES} )
#SET_TARGET_PROPERTIES ( GBTL_algorithms_static PROPERTIES VERSION 1.2.0
#    OUTPUT_NAME graphblas
#    POSITION_INDEPENDENT_CODE OFF
#    SOVERSION 1
#    C_STANDARD_REQUIRED 11
#    PUBLIC_HEADER "GBTL/algorithms.h" )
#set_property ( TARGET graphblas_static PROPERTY C_STANDARD 11 )

# graphblas installation location
#install ( TARGETS graphblas graphblas_static
#    LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#    ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} )

# Demo library
#file ( GLOB DEMO_SOURCES "Demo/Source/*.c" )
#add_library ( graphblasdemo SHARED ${DEMO_SOURCES} )
#SET_TARGET_PROPERTIES ( graphblasdemo PROPERTIES
#    C_STANDARD_REQUIRED 11 )
#set_property ( TARGET graphblasdemo PROPERTY C_STANDARD 11 )

#add_library ( graphblasdemo_static STATIC ${DEMO_SOURCES} )
#SET_TARGET_PROPERTIES ( graphblasdemo_static PROPERTIES
#    C_STANDARD_REQUIRED 11 )
#set_property ( TARGET graphblasdemo_static PROPERTY C_STANDARD 11 )

# Demo programs
add_executable ( hpec_ex0 "tutorial/hpec_ex0.c" )
add_executable ( hpec_ex1_hello_world "tutorial/hpec_ex1_hello_world.c" )
add_executable ( hpec_ex2_neighbors   "tutorial/hpec_ex2_neighbors.c" )
add_executable ( hpec_ex3_level_bfs   "tutorial/hpec_ex3_level_bfs.c" )
add_executable ( hpec_ex4_bfs_parents "tutorial/hpec_ex4_bfs_parents.c" )

# Libraries required for Demo programs
target_link_libraries ( hpec_ex0             graphblas )
target_link_libraries ( hpec_ex1_hello_world graphblas )
target_link_libraries ( hpec_ex2_neighbors   graphblas )
target_link_libraries ( hpec_ex3_level_bfs   graphblas )
target_link_libraries ( hpec_ex4_bfs_parents graphblas )

