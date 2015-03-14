// -----------------------------------------------------------------------
//
//                                gaussian_fit4_nlls.cpp V 0.01
//
//                                 (c) Brian Lynch March, 2015
//
// -----------------------------------------------------------------------

#include <math.h>
#include <vector>
#include <iostream>
#include <new>

#include "gaussian_fit4_nlls.h"
#include "lif_analysis.h"
#include "matrix_utils/matrix_ops.h"

/************************************************************************/
/*
 * 4 parameter nonlinear least squares fitting:
 *      http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html
 */
int gauss_fit4_nlls(double **x, double **fx,
                    const unsigned int &Npoints, const unsigned int &Ntries,
                      const double &TOL, struct GaussFit4Params &FitParams){
//std::cout << "BEGIN gaussian_fit4_nlls" << std::endl;

   int res  = 0;
   
   unsigned int it   = 0,
                Npar = FitParams.Npar; //# fit parameters [4 in this case]
       
   double xt      = 0.0,  // Temporary Counts
          fxt     = 0.0,  // Temporary Lambda
          *dFx    = NULL, // Difference between fit and data
          *R2     = NULL, // Sum of squared residuals
          *dparam = NULL, // Difference between new and old parameters
          *A      = NULL, // A matrix
          *AT     = NULL, // Tranposed A matrix
          *a      = NULL, // Product of AT * A
          *ainv   = NULL, // Inverse of AT * A
          *I      = NULL, // Identity matrix can be used for diagnostics
          *b      = NULL, // Product of AT * dFx
          *param  = NULL; // Parameter array storing struc IVFIT2Params info
      
   // Function pointer to help setup the A and AT matrices
   double (*FXds[])(const double &, const double &, const double &,
                                  const double &, const double &) =
                                   {dFxdxo, dFxdsig2, dFxdA, dFxdB};
   
   // Since there may be ALOT of data points, make sure new[] is successful
   try{
      
      dFx    = new double[Npoints],
      R2     = new double[Ntries],
      dparam = new double[Npar],
      A      = new double[Npoints * Npar],
      AT     = new double[Npar * Npoints],
      a      = new double[Npar * Npar],
      ainv   = new double[Npar * Npar],
      I      = new double[Npar * Npar],
      b      = new double[Npar],
      param  = new double[Npar];
      
   }catch(std::bad_alloc& ba){
      
      std::cerr << "ERROR: gaussian_fit4_nlls initialization: " << ba.what();
      std::cerr << std::endl;
      res = 0;
      goto cleanup;
      
   }
   
   //Set the initial fit parameters
   param[0] = FitParams.x0;
   param[1] = FitParams.sigma2;
   param[2] = FitParams.Ao;
   param[3] = FitParams.Bo;
   R2[0] = 1.0;
    
   while((it < Ntries) && (R2[it] > TOL)){
      
      //Calculate the A Matrix
      for(unsigned int row = 0; row < Npoints; row++){
         
         xt       = (*x)[row];
         fxt      = Fxa(xt,param[0],param[1],param[2],param[3]);
         dFx[row] = (*fx)[row] - fxt;
         
         for(unsigned int col = 0; col < Npar; col++){
            
            A[row * Npar + col] = FXds[col](xt,param[0],param[1],param[2],param[3]);
            
         }

      }
      //getchar();
      
      //Now find the transpose of A
      if(!TransposeMatrix(A, Npoints, Npar, &AT)){
       
         std::cerr << "ERROR: transposing matrix failed: A" << std::endl;
         res = 0;
         goto cleanup;
         
      }
       
      //Product of a = AT * A. a is the matrix we need to invert
      if(!MultiplyMatrix(AT, Npar, Npoints, A, Npoints, Npar, &a)){
       
         std::cerr << "ERROR: matrix multiplication failed: AT * A";
         std::cerr << std::endl;
         res = 0;
         goto cleanup;
         
      }
      
      //Calculate the inverse matrix ainv
      if(!InvertMatrix(a, Npar, &ainv)){
       
         std::cerr << "ERROR: matrix inversion failed: ainv" << std::endl;
         PrintMatrix(ainv,Npar,Npar);
         std::cerr << "\n";
         PrintMatrix(a,Npar,Npar);
         res = 0;
         goto cleanup;
         
      }
      
      /*
       * Calculate the product of ainv * a [should be the identity matrix].
       * If trouble with routine, could be used as diagnostic.
       */
      if(!MultiplyMatrix(ainv, Npar, Npar, a, Npar, Npar, &I)){
       
         std::cerr << "ERROR: matrix multiplication failed: ainv * a";
         std::cerr << std::endl;
         PrintMatrix(I, Npar, Npar);
         res = 0;
         goto cleanup;
         
      }
      
      /*
       * Product of AT and the difference between data and model.
       */
      if(!MultiplyMatrix(AT, Npar, Npoints, dFx, Npoints, 1, &b)){
       
         std::cerr << "ERROR: matrix multiplication failed: AT * dIi";
         std::cerr << std::endl;
         res = 0;
         goto cleanup;
         
      }
    
      //Calculate the small increment toward convergence
      if(!MultiplyMatrix(ainv, Npar, Npar, b, Npar, 1, &dparam)){
       
         std::cerr << "ERROR: matrix multiplication failed:  ainv * b";
         std::cerr << std::endl;
         res = 0;
         goto cleanup;
         
      }
    
      /* 
       * Final update tasks:
       *        1) Increment the iteration #
       *        2) New param values by applying offset
       *        3) The sum of squared residuals to check convergence
       */
      ++it;
      R2[it] = 0.0;
      for(unsigned int i = 0; i < Npar; i++){
         
         param[i] += dparam[i];
         R2[it] += dparam[i] * dparam[i];  
         
      }
  
   }//End while loop checking convergence tolerance or max iterations
   
   res = 1;
   
   //Print Results
   //Set the initial fit parameters
   FitParams.x0     = param[0];
   FitParams.sigma2 = param[1]; 
   FitParams.Ao     = param[2];
   FitParams.Bo     = param[3];
   std::cout << " R^2         : " << R2[it] << std::endl;
   std::cout << " # iterations: " << it << std::endl;
   
//Memory cleanup
cleanup:
   
   delete[] dFx;
   delete[] R2;
   delete[] dparam;
   delete[] A;
   delete[] AT;
   delete[] a;
   delete[] ainv;
   delete[] I;
   delete[] b;
   delete[] param;
  
//std::cout << "END gaussian_fit4_nlls" << std::endl;
return (res);
}//End function gaussian_fit4_nlls