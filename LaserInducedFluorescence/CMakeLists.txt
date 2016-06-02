# ------------------------------------------------------------------------
#
#                            CMakeLists.txt for LIFAnalysis
#                                        V 0.01
#
#                             (c) Brian Lynch March, 2015
#
# ------------------------------------------------------------------------
cmake_minimum_required (VERSION 2.6)
project(LIFAnalysis)

#Tell cmake to look in the following subdirectories
#for other files named CMakeLists.txt
add_subdirectory (lif)
add_subdirectory (matrix_utils)