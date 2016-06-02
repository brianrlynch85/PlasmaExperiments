# ------------------------------------------------------------------------
#
#                         CMakeLists.txt for the matrix_utils
#                                        V 0.01
#
#                            (c) Brian Lynch February, 2015
#
# ------------------------------------------------------------------------
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib)
add_library(matrix_utilslib matrix_ops.cpp)