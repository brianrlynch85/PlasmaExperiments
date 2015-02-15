// -----------------------------------------------------------------------------------------------
//
//                                     IVFit2NLLS.h V 0.01
//
//                                (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------------------------------

#ifndef IVFit2NLLS_h
#define IVFit2NLLS_h

#include "DoubleProbeAnalysis.h"

/************************************************************************/
/*
 * IVFIT2NLLS(...) performs a 2 parameter non-linear least squares curve fit:
 *       
 *      @param[in] std::vector Ii: an input vector of current measurements
 *      @param[in] std::vector V: an input vector of voltage measurements
 *      @param[in] int Ntries: maximum # attempts to curve fit 
 *      @param[in] double TOLERANCE: convergence tolerance
 *      @param[in/out] struct IVFit2Params Fitparams: input guess / output final fit paramters
 *      @return int success/failure
 * 
 */
int IVFit2NLLS(const std::vector<double> &Ii, const std::vector<double> &V,
                       const unsigned int &Ntries, const double &TOLERANCE,
                                           struct IVFit2Params &FitParams);

/************************************************************************/
/*
 * The typical double probe characteristic trace is given by:
 * 
 *      I(V) = Isat * tanh(0.5 * e * V / Te)
 *              
 *      @param[in] double V: the voltage difference between the probe tips
 *      @param[in] double e: the electron charge
 *      @param[in] double Te: the electron temperature
 *      @return double: the collected current
 * 
 */
inline double Iv(const double &V, const double &Isat, const double &Te){
   
   double result = 0.0;
   
   result = Isat * tanh(0.5 * V / Te);
      
   return (result);
      
}

/************************************************************************/
/*
 * The partial derivative of I(V) w.r.t Isat:
 * 
 *      d(I(V))/d(Isat) = tanh(0.5 * e * V / Te)
 *              
 *      @param[in]  V: the voltage difference between the probe tips
 *      @param[in]  e: the electron charge
 *      @param[in] Te: the electron temperature
 *      @return double: the collected current
 * 
 */
inline double dIvdIsat(const double &V, const double &Isat, const double &Te){
  
   double result = 0.0;
   
   result = tanh(0.5 * V / Te);
   
   return (result);
   
}

/************************************************************************/
/*
 * The partial derivative of I(V) w.r.t Te:
 * 
 *      d(I(V))/d(Te) = - Isat * 0.5 * V * (1.0 - pow(tanh(0.5 * V / Te),2.0)) / pow(Te,2.0)
 *              
 *      @param[in]  V: the voltage difference between the probe tips
 *      @param[in]  e: the electron charge
 *      @param[in] Te: the electron temperature
 *      @return double: the collected current
 * 
 */
inline double dIvdTe(const double &V, const double &Isat, const double &Te){
  
   double result = 0.0;
   
   result  = Isat * 0.5 * V * (1.0 - pow(tanh(0.5 * V / Te),2.0));
   result *= -1.0 / (Te * Te);
   
   return (result);
   
}

#endif
