// -----------------------------------------------------------------------------------------------
//
//                                     IVFit2NLLS.cpp V 0.01
//
//                                 (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------------------------------

#include <math.h>
#include <vector>
#include <iostream>
#include <new>

#include "IVFit2NLLS.h"
#include "DoubleProbeAnalysis.h"
#include "matrix_utils/matrix_ops.h"

/************************************************************************/
/*
 * 2 parameter nonlinear least squares fitting:
 *      http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html
 */
int IVFit2NLLS(const std::vector<double> &Ii, const std::vector<double> &V,
                       const unsigned int &Ntries, const double &TOLERANCE,
                                           struct IVFit2Params &FitParams){
//std::cout << "BEGIN IVFit2NLLS" << std::endl;

   int res  = 0;
   
   unsigned int it   = 0,
                Npar = FitParams.Npar, //# fit parameters [2 in this case]
                Npoi = V.size();       //# data points stored in I and V
       
   double Ivt     = 0.0,  //Temporary collected current
          Vt      = 0.0,  //Temporary voltage measurement
          *dIi    = NULL, //Difference between fit and data
          *R2     = NULL, //Sum of squared residuals
          *dparam = NULL, //Difference between new and old parameters
          *J      = NULL, //Jacobian matrix
          *JT     = NULL, //Tranposed Jacobian matrix
          *a      = NULL, //Product of JT * J
          *ainv   = NULL, //Inverse of JT * J
          *I      = NULL, //Identity matrix can be used for diagnostics
          *b      = NULL, //Product of JT * dIi
          *param  = NULL; //Parameter array storing struc IVFIT2Params info
      
   //Function pointer to help setup the Jacobian
   double (*IVds[])(const double &,const double &, const double &) = {dIvdIsat, dIvdTe};
   
   //Make sure number of read points for Ii and V are the same before proceeding
   if(Ii.size() != V.size()){
      
      std::cout << "Passed incompatible arrays for I and V input data" << std::endl;
      res = 0;
      return(res);
     
   }
   
   //Since there could be ALOT of data points, make sure new[] is successful
   try{
      
      dIi    = new double[Npoi],
      R2     = new double[Ntries],
      dparam = new double[Npar],
      J      = new double[Npoi * Npar],
      JT     = new double[Npar * Npoi],
      a      = new double[Npar * Npar],
      ainv   = new double[Npar * Npar],
      I      = new double[Npar * Npar],
      b      = new double[Npar],
      param  = new double[Npar];
      
   }catch(std::bad_alloc& ba){
      
      std::cerr << "ERROR: in IVFIT2NLLS initialization: " << ba.what() << std::endl;
      res = 0;
      goto cleanup;
      
   }
   
   //Set the initial fit parameters
   param[0] = FitParams.Isat;
   param[1] = FitParams.Te;
   R2[0] = 1.0;
    
   while((it < Ntries) && (R2[it] > TOLERANCE)){
   
      for(unsigned int row = 0; row < Npoi; row++){
         
         Vt = V.at(row);
         Ivt = Iv(Vt,param[0],param[1]);
         dIi[row] = (Ii.at(row)-Ivt);
         
         for(unsigned int col = 0; col < Npar; col++){
            
            //Calculate the Jacobian Matrix
            J[row * Npar + col] = IVds[col](Vt,param[0],param[1]);
            
         }

      }
      
      for(unsigned int row = 0; row < Npar; row++){
         
         for(unsigned int col = 0; col < Npoi; col++){
            
            Vt = V.at(col);
            
            //Calculate the transpose Jacobian Matrix
            JT[row * Npoi + col] = IVds[row](Vt,param[0],param[1]);
            
         }

      }
      
      /* 
       * Product of transposed Jacobian with Jacobian. This
       * is the matrix we need to invert.
       */
      if(!MultiplyMatrix(JT, Npar, Npoi, J, Npoi, Npar, &ainv)){
       
         std::cout << "ERROR: matrix multiplication failed" << std::endl;
         res = 0;
         goto cleanup;
         
      }
      
      //Just have have a copy of the uninverted matrix
      if(!MultiplyMatrix(JT, Npar, Npoi, J, Npoi, Npar, &a)){
       
         std::cout << "ERROR: matrix multiplication failed" << std::endl;
         res = 0;
         goto cleanup;
         
      }
      
      //Calculate the inverse matrix ainv
      if(!InvertMatrix(&ainv, Npar)){
       
         std::cout << "ERROR: matrix inversion failed" << std::endl;
         res = 0;
         goto cleanup;
         
      }
      
      /*
       * Calculate the product of ainv * a [result should be identity matrix].
       * If trouble with routine, could be used as diagnostic.
       */
      if(!MultiplyMatrix(ainv, Npar, Npar, a, Npar, Npar, &I)){
       
         std::cout << "ERROR: matrix multiplication failed" << std::endl;
         res = 0;
         goto cleanup;
         
      }
      
      /*
       * Product of tranposed Jacobian and the difference between
       * data and model.
       */
      if(!MultiplyMatrix(JT, Npar, Npoi, dIi, Npoi, 1, &b)){
       
         std::cout << "ERROR: matrix multiplication failed" << std::endl;
         res = 0;
         goto cleanup;
         
      }
    
      //Calculate the small increment toward convergence
      if(!MultiplyMatrix(ainv, Npar, Npar, b, Npar, 1, &dparam)){
       
         std::cout << "ERROR: matrix multiplication failed" << std::endl;
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
   FitParams.Isat = param[0];
   FitParams.Te   = param[1];
   std::cout << " R^2         : " << R2[it] << std::endl;
   std::cout << " # iterations: " << it << std::endl;
   
//Memory cleanup
cleanup:
   
   delete[] dIi;
   delete[] R2;
   delete[] dparam;
   delete[] J;
   delete[] JT;
   delete[] a;
   delete[] ainv;
   delete[] I;
   delete[] b;
   delete[] param;
  
//std::cout << "END IVFit2NLLS" << std::endl;
return (res);
}//End function IVFit2NNLS