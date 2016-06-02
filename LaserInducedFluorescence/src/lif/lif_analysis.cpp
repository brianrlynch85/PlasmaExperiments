// -----------------------------------------------------------------------
//
//                                   lif_analysis.cpp V 0.01
//
//                                 (c) Brian Lynch March, 2015
//
// -----------------------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <getopt.h>
#include <iomanip>

#include "gaussian_fit4_nlls.h"
#include "lif_analysis.h"

/************************************************************************/
int main(int argc, char** argv){

std::cout << "-- BEGIN lif_analysis --" << std::endl;

   /*
    * Initial guess parameters for the Non-linear least squares
    * fit (4 parameters).
    */
   const double xo_guess   = 668.6138;  // Rest wavelength [nm]
   const double sig2_guess = 0.0000006; // Standard deviation [nm^2]
   const double Ao_guess   = 4.0;       // Signal amplitude
   const double Bo_guess   = 0.5;       // Signal background

   const int    Max = 100;    // Maximum number of iterations while fitting
   const double Tol = 1.0E-8; // Tolerance for convergence of the curve fit

   int opt = 0;                 // Command line option parser variable
   char *input_filename = NULL; // Command line option input file

   // Parse the command line
   if(1 == argc){
      
      print_usage(); 
      exit(EXIT_FAILURE);
       
   }
      
   while((opt = getopt(argc, argv,"-f:")) != -1) {
     
      switch (opt) {
         
         case 'f' : // Input filename option
            
            input_filename = optarg;
            std::cout << "Input Filename: " << input_filename;
            std::cout << std::endl;
            break;
            
         case '?': // Unrecognized command line option
            
            std::cerr << "Unrecognized command line option";
            std::cerr << std::endl;
            print_usage();
            return (-1);
            
         case '\1': // The - in "-f:" finds non option command line params
            
            std::cerr << "Passed non-option to command line";
            std::cerr << std::endl;
            print_usage();
            return (-1);
            
      }
        
   }
   
   std::vector<double> lambda, // Input wavelength
                       counts; // Input counts
                       
   double *la = NULL, // Array for input lambdas
          *ca = NULL; // Array for input counts
          
   unsigned int Na = 0; // Size of input arrays

   // Array used to store initial fit parameter guesses
   struct GaussFit4Params FitParams = {xo_guess, sig2_guess, Ao_guess, Bo_guess, 4};
   std::cout.precision(7);
   std::cout << "Initial fit parameters: " << std::endl;
   std::cout << " Rest Wavelength        [nm]  : " << xo_guess << std::endl;
   std::cout << " Sigma^2               [nm^2] : " << sig2_guess << std::endl;
   std::cout << " Amplitude               []   : " << Ao_guess << std::endl;
   std::cout << " Background              []   : " << Bo_guess << std::endl;
   
   // Input and Output string information
   std::string input_filename_s(input_filename);
   std::ifstream input_file(input_filename_s.c_str(), std::ifstream::in);
   
   // Attempt to read the input file
   if(input_file.is_open()){
     
      std::cout << "Reading data..." << std::endl;
      double col1 = 0.0,
             col2 = 0.0;
     
      while(input_file.good()){
        
         input_file >> col1 >> col2;
        
         if(!input_file.eof()){ //Last line is not stored twice
           
            lambda.push_back(col1);
            counts.push_back(col2);
           
         }
       
      }
     
   }else{
     
      std::cerr << "Error opening file:" << input_filename_s.c_str();
      std::cerr << std::endl;
      return (-1);
      
   }//Done attempting to read file
   
   input_file.close();
   
   //Copy the vectors used to read file into arrays
   Na = lambda.size();
   la = new double[Na];
   ca = new double[Na];
   std::copy(lambda.begin(), lambda.end(), la);
   std::copy(counts.begin(), counts.end(), ca);
      
   //Perform the double probe curve fit using non-linear least squares
   std::cout << "Performing curve fit..." << std::endl;
   if(gauss_fit4_nlls(&la, &ca, Na, Max, Tol, FitParams)){
      
      std::cout << "Curve fit successful!" << std::endl;
      
   }else{
    
      std::cout << "Curve fit FAILED!" << std::endl;
      
   }
   
   //Print the fitted parameters
   std::cout.precision(7);
   std::cout << "Final fit parameters: " << std::endl;
   std::cout << " Rest Wavelength        [nm]  : " << FitParams.x0 << std::endl;
   std::cout << " Sigma^2               [nm^2] : " << FitParams.sigma2 << std::endl;
   std::cout << " Amplitude               []   : " << FitParams.Ao << std::endl;
   std::cout << " Background              []   : " << FitParams.Bo << std::endl;
   
   //Declare and write the output file
   std::string output_filename_s(input_filename);
   output_filename_s.resize(output_filename_s.length()-4);
   output_filename_s.append("_fit.dat");
   std::ofstream output_file(output_filename_s.c_str(), std::ofstream::out);
   output_file << std::scientific;
   
   //Attempt to open the output file
   if(output_file.is_open()){
     
      std::cout << "Writing fit data to file: " << output_filename_s.c_str();
      std::cout << std::endl;
      double col1       = lambda.at(0),
             col2       = 0.0,
             lambda_end = *std::max_element(lambda.begin(),lambda.end());
     
      while(col1 < lambda_end){
            
         // Since the input data scans back and forth, lets only
         // use the values from forward back of scan.
         col2 = Fxa(col1, FitParams.x0, FitParams.sigma2, FitParams.Ao,
                                                            FitParams.Bo);
         output_file << col1 << " ";
         output_file << col2 << std::endl;
         col1 += 0.0001;
         
      }
     
   }else{
     
      std::cerr << "Error opening file:" << output_filename_s.c_str();
      std::cerr << std::endl;
      return (-1);
      
   }//Done writing output file
   
   output_file.close();

   delete[] ca;
   delete[] la;
 
std::cout << "-- END lif_analysis --" << std::endl;
return(0);

}

/************************************************************************/
void print_usage(){
   
   std::cout << "Usage:" << std::endl;
   std::cout << "build/bin/LIFAnalysis -f <filename>" << std::endl;
   std::cout << "build/bin/LIFAnalysis -f ExampleData/ExampleData.dat";
   std::cout << std::endl;
   
}
