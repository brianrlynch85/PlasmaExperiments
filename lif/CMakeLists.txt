# ------------------------------------------------------------------------
#
#                                 CMakeLists.txt for lif
#                                        V 0.01
#
#                            (c) Brian Lynch March, 2015
#
# ------------------------------------------------------------------------
cmake_minimum_required (VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)

#Set gdb and warning flags
set(CMAKE_CXX_FLAGS "-Wall")

#Make sure lapack is installed
find_package(LAPACK REQUIRED)

#Incluce this directory
include_directories(${PROJECT_SOURCE_DIR})

#Print the included directories
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("inc_dirs = ${inc_dirs}")

#Set the executable lif_analysis source dependencies
set(lif_src lif_analysis.cpp gaussian_fit4_nlls.cpp)

#Add the executable, which will be in build/bin
add_executable(LIFAnalysis ${lif_src})

#Link the matrix_utils library and LAPACK
target_link_libraries(LIFAnalysis matrix_utilslib)
target_link_libraries(LIFAnalysis ${LAPACK_LIBRARIES})