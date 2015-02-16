// -----------------------------------------------------------------------
//
//                                    matrix_ops.cpp V 0.01
//
//                                (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "matrix_ops.h"

/************************************************************************/
/*
 * The functions in this header may make use of the  LAPACK Linear Algebra
 * library. For more information, see LAPACK: http://www.netlib.org/lapack/
 * as well as the license file llapack_license contained in the
 * matrix_utils folder.
 * 
 * A is a N x N matrix stored as a linear data structure
 * 
*/
int InvertMatrix(double **A, int N){
   
   int res = 0,
       INFO   = -1,        //Status helper from lapack functions
       LWORK = N * N * N;  //How big we need the work space to be
       
   //Pivot indices of the matrix for row swaps used in inversion
   int *IPIV = (int *)malloc((N+1) * sizeof(int));
    
   //Workspace declaration from Lwork
   double *WORK = (double *)malloc(LWORK * sizeof(double));
   
   if(NULL == IPIV){
      
      printf("malloc of *IPIV pointer failed inside InvertMatrix\n");
      res = 0;
      goto cleanup;
      
   }
   
   if(NULL == WORK){
      
      printf("malloc of *WORK pointer failed inside InvertMatrix\n");
      res = 0;
      goto cleanup;
      
   }

   // LU decomoposition of a general matrix
   //http://www.netlib.no/netlib/lapack/double/dgetrf.f
   dgetrf_(&N,&N,*A,&N,IPIV,&INFO);
   
   if(INFO != 0){
       
      printf("Problem with LU decompositio inside function dgetrf\n");
      printf("(see LAPACK documentation): INFO = %d\n",INFO);  
      res= 0;
      goto cleanup;
       
   }else{
      
      //printf("Matrix LU factorization successful\n");
      
   }
   
   //Find the inverse of a matrix A given its LU decomposition
   //http://www.netlib.no/netlib/lapack/double/dgetri.f
   dgetri_(&N,*A,&N,IPIV,WORK,&LWORK,&INFO);
   
   if(INFO != 0){
       
      printf("Problem with matrix inversion inside function dgetri\n");
      printf("(see LAPACK documentation): INFO = %d\n",INFO);  
      res = 0;
      goto cleanup;
       
   }else{
      
      res = 1;
      //printf("Matrix inversion successful\n");
      
   }

//Memory cleanup
cleanup:
   
   free(IPIV);
   free(WORK);
   
return(res);
}; //End function InvertMatrix

/************************************************************************/
/* 
 * Function multiplies 2 matrices C = A * B
 * It is the users responsibility to make sure C is properly allocated.
 */
int MultiplyMatrix(const double *A, const int &ANROW, const int &ANCOL,
                   const double *B, const int &BNROW, const int &BNCOL,
                                                           double **C){
   
   int res = 0;
   
   double dot = 0.0;
   
   //Make sure the matrices can actually be multiplied
   if(ANCOL != BNROW){
    
      printf("ERROR: MultiplyMatrix requires # A Cols = # B Rows\n");
      res = 0;
      return (res);
      
   }
   
   for(int row = 0; row < ANROW; row++){
         
      for(int col = 0; col < BNCOL; col++){
            
         dot = 0.0;
            
         for(int l = 0; l < ANCOL; l++){
              
            dot += A[row * ANCOL + l] * B[l * BNCOL + col];
              
         }
         
         (*C)[row * BNCOL + col] = dot;

      }
         
   }
   
   res = 1;
   
return(res);
} //End function MultiplyMatrix