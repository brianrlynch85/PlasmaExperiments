# ------------------------------------------------------------------------
#
#                            CMakeLists.txt for doubleprobe
#                                        V 0.01
#
#                            (c) Brian Lynch February, 2015
#
# ------------------------------------------------------------------------
cmake_minimum_required (VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)

#Set gdb and warning flags
set(CMAKE_CXX_FLAGS "-g -Wall")

#Make sure lapack is installed
find_package(LAPACK REQUIRED)

#Incluce this directory
include_directories(${CMAKE_SOURCE_DIR})

#Print the included directories
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("inc_dirs = ${inc_dirs}")

#Set the executable DoubleProbeAnalysis source dependencies
set(dpa_src DoubleProbeAnalysis.cpp IVFit2NLLS.cpp)

#Add the executable, which will be in build/bin
add_executable(DoubleProbeAnalysis ${dpa_src})

#Link the matrix_utils library and LAPACK
target_link_libraries(DoubleProbeAnalysis matrix_utilslib)
target_link_libraries(DoubleProbeAnalysis ${LAPACK_LIBRARIES})