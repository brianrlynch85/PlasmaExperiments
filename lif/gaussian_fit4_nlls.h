// -----------------------------------------------------------------------
//
//                                gaussian_fit4_nlls.h V 0.01
//
//                                (c) Brian Lynch March, 2015
//
// -----------------------------------------------------------------------

#ifndef lif_gaussian_fit4_nlls_h
#define lif_gaussian_fit4_nlls_h

#include <stdio.h>
#include "lif_analysis.h"

/************************************************************************/
/*
 * gauss_fit4_nlls(...) performs a 4 parameter nonlinear least squares
 * curve fit:
 *       
 *      @param[in] x            : input array of wavelengths
 *      @param[in] fx           : input array of # counts
 *      @param[in] Npoints      : length of input arrays
 *      @param[in] Ntries       : maximum # attempts to curve fit 
 *      @param[in] TOL          : convergence tolerance
 *      @param[in/out] Fitparams: input guess / output final fit paramters
 *      @return int success/failure
 * 
 */
int gauss_fit4_nlls(double **x, double **fx,
                    const unsigned int &Npoints, const unsigned int &Ntries,
                      const double &TOL, struct GaussFit4Params &FitParams);

/************************************************************************/
/*
 * The typical LIF characteristic trace is given by:
 * 
 *      = A * exp(-0.5 * (x - xo) * (x - xo) / sig2) + B
 *     
 *      @param[in] x    : independent variable
 *      @param[in] xo   : mean
 *      @param[in] sig2 : standard deviation
 *      @param[in] A    : amplitude
 *      @param[in] B    : background noise offset
 *      @return Fxa
 * 
 */
inline double Fxa(const double &x, const double &xo, const double &sig2,
                                      const double &A, const double &B){
       
   return (A * exp(-0.5 * (x - xo) * (x - xo) / sig2) + B);
      
}

/************************************************************************/
/*
 * Partial derivative of Fxa w.r.t. Ao:
 * 
 *      = exp(-0.5 * (x - x0) * (x - x0) / sig2)
 *     
 *      @param[in] x    : independent variable
 *      @param[in] xo   : mean
 *      @param[in] sig2 : standard deviation
 *      @param[in] A    : amplitude
 *      @param[in] B    : background noise offset
 *      @return dFxdA 
 * 
 */
inline double dFxdA(const double &x, const double &xo, const double &sig2,
                                        const double &A, const double &B){
   
   return (exp(-0.5 * (x - xo) * (x - xo) / sig2));
   
}

/************************************************************************/
/*
 * Partial derivative of Fxa w.r.t. xo:
 * 
 *      = A * (x - xo) * exp(-0.5 * (x - xo) * (x - xo) / sig2) / sig2
 *     
 *      @param[in] x    : independent variable
 *      @param[in] xo   : mean
 *      @param[in] sig2 : standard deviation
 *      @param[in] A    : amplitude
 *      @param[in] B    : background noise offset
 *      @return dFxdxo 
 * 
 */
inline double dFxdxo(const double &x, const double &xo, const double &sig2,
                                         const double &A, const double &B){
   
   return (A * (x - xo) * exp(-0.5 * (x - xo) * (x - xo) / sig2) / sig2);
   
}

/************************************************************************/
/*
 * Partial derivative of Fxa w.r.t. sig2:
 * 
 *      = A * 0.5 * (x - xo) * (x - xo) * exp(-0.5 * (x - xo) * (x - xo) / sig2) / (sig2 * sig2)
 *     
 *      @param[in] x    : independent variable
 *      @param[in] xo   : mean
 *      @param[in] sig2 : standard deviation
 *      @param[in] A    : amplitude
 *      @param[in] B    : background noise offset
 *      @return dFxdsig2
 * 
 */
inline double dFxdsig2(const double &x, const double &xo, const double &sig2,
                                           const double &A, const double &B){
   
   return (A * 0.5 * (x - xo) * (x - xo) * exp(-0.5 * (x - xo) * (x - xo) / sig2) / (sig2 * sig2));
   
}

/************************************************************************/
/*
 * Partial derivative of Fxa w.r.t. Bo:
 * 
 *      = 1.0
 *     
 *      @param[in] x    : independent variable
 *      @param[in] xo   : mean
 *      @param[in] sig2 : standard deviation
 *      @param[in] A    : amplitude
 *      @param[in] B    : background noise offset
 *      @return dFxdB
 * 
 */
inline double dFxdB(const double &x, const double &xo, const double &sig2,
                                        const double &A, const double &B){
   
   return (1.0);
   
}

#endif