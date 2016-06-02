// -----------------------------------------------------------------------
//
//                                  DoubleProbeAnalysis.h V 0.01
//
//                                 (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------

#ifndef DoubleProbeAnalysis_h
#define DoubleProbeAnalysis_h

struct IVFit2Params{
  
   double Isat; //Ion saturation current [A]
   double Te;   //Electron temperature [eV]
   int    Npar; //Number of parameters (2)
   
};

/************************************************************************/
/*
 * Usage function used to display example calling commands.
 */
void print_usage();

#endif