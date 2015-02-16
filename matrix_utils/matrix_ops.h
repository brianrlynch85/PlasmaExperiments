// -----------------------------------------------------------------------
//
//                                     matrix_ops.h V 0.01
//
//                                (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------

#ifndef matrix_ops_h
#define matrix_ops_h

/************************************************************************/
/*
 * The functions in this header may make use of the  LAPACK Linear Algebra
 * library. For more information, see LAPACK: http://www.netlib.org/lapack/
 * as well as the license file llapack_license contained in the
 * matrix_utils folder.
 * 
 * Wrap the fortran functions in extern so a C++ compiler does not load
 * argument information ("name mangling") when calling fortran functions
 * dgetrf and dgetri.
 * 
*/
extern "C"{
   
    //Find the LU decomoposition of matrix A with dimension M x N
   //http://www.netlib.no/netlib/lapack/double/dgetrf.f
    void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV,
                                                     int *INFO);

    //Find the inverse of a matrix A given its LU decomposition
    //http://www.netlib.no/netlib/lapack/double/dgetri.f
    void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK,
                                               int *lwork, int *INFO);
    
}

/************************************************************************/
/*
 * InvertMatrix(...) calculates the inverse of A. A is a N x N matrix
 * stored as a linear data structure. 
 *              
 *      @param[in/out] double **A: matrix A
 *      @param[in] int N: # rows and cols in A
 *      @return int: success/failure
 * 
 */
int InvertMatrix(double **A, int N);

/************************************************************************/
/*
 * MultiplyMatrix(...) calculates the product A * B = C
 *              
 *      @param[in] double *A: matrix A
 *      @param[in] int ANROW: # rows in A
 *      @param[in] int ANCOL: # cols in A
 *      @param[in] double *B: matrix B
 *      @param[in] int BNROW: # rows in B
 *      @param[in] int BNCOL: # cols in B
 *      @param[in] double **C: matrix C
 *      @return int: success/failure
 * 
 */
int MultiplyMatrix(const double *A, const int &ANROW, const int &ANCOL,
                   const double *B, const int &BNROW, const int &BNCOL,
                                                           double **C);

#endif
