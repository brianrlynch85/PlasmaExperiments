// -----------------------------------------------------------------------
//
//                                     lif_analysis.h V 0.01
//
//                                 (c) Brian Lynch March, 2015
//
// -----------------------------------------------------------------------

#ifndef lif_lif_analysis_h
#define lif_lif_analysis_h

struct GaussFit4Params{
  
   double x0     ; // Rest wavelength              [m]
   double sigma2 ; // Ion temperature              [eV]
   double Ao     ; // Amplitude of arbitary counts []
   double Bo     ; // Amplitude of background      []
   int    Npar   ; // Number of parameters (4)
   
};

/************************************************************************/
/*
 * Usage function used to display example calling commands.
 */
void print_usage();

#endif