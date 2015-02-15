// -----------------------------------------------------------------------------------------------
//
//                                 DoubleProbeAnalysis.cpp V 0.01
//
//                                 (c) Brian Lynch February, 2015
//
// -----------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <getopt.h>
#include <iomanip>

#include "IVFit2NLLS.h"
#include "DoubleProbeAnalysis.h"

/************************************************************************/
int main(int argc, char** argv){

std::cout << "-- BEGIN DoubleProbeAnalysis --" << std::endl;

   /*
    * Initial guess parameters for the Non-linear least squares
    * fit (2 parameters).
    */
   const double Is_guess = 3.3E-6; //Ion saturation current [A]
   const double Te_guess = 3.0;    //Electron temperature   [eV]

   const int    Max = 100;    //Maximum number of iterations while fitting
   const double Tol = 1.0E-8; //Tolerance for convergence of the curve fit

   int opt = 0;                 //Command line option parser variable
   char *input_filename = NULL; //Command line option input file

   //Parse the command line
   if(1 == argc){
      
      print_usage(); 
      exit(EXIT_FAILURE);
       
   }
      
   while((opt = getopt(argc, argv,"-f:")) != -1) {
     
      switch (opt) {
         
         case 'f' : //input filename option
            
            input_filename = optarg;
            std::cout << "Input Filename: " << input_filename;
            std::cout << std::endl;
            break;
            
         case '?': //unrecognized command line option
            
            std::cerr << "Unrecognized command line option";
            std::cerr << std::endl;
            print_usage();
            return (-1);
            
         case '\1': //the - in "-f:" finds non option command line params
            
            std::cerr << "Passed non-option to command line";
            std::cerr << std::endl;
            print_usage();
            return (-1);
            
      }
        
   }
   
   std::vector<double> Ii, //I (input current)
                       Vi, //V (input voltage)
                       Ia; //I (current using fitted parameters)

   //Array used to store initial fit parameter guesses
   struct IVFit2Params FitParams = {Is_guess, Te_guess, 2};
   std::cout.precision(3);
   std::cout << "Initial fit parameters: " << std::endl;
   std::cout << " Ion saturation current [A]  : " << Is_guess << std::endl;
   std::cout << " Electron temperature   [eV] : " << Te_guess << std::endl;
   
   //Input and Output string information
   std::string input_filename_s(input_filename);
   std::ifstream input_file(input_filename_s.c_str(), std::ifstream::in);
   
   //Attempt to read the input file
   if(input_file.is_open()){
     
      std::cout << "Reading IV data..." << std::endl;
      double col1 = 0.0,
             col2 = 0.0;
     
      while(input_file.good()){
        
         input_file >> col1 >> col2;
        
         if(!input_file.eof()){ //Last line is not stored twice
           
            Vi.push_back(col1);
            Ii.push_back(col2);
           
         }
       
      }
     
   }else{
     
      std::cerr << "Error opening file:" << input_filename_s.c_str();
      std::cerr << std::endl;
      return (-1);
      
   }//Done attempting to read file
   
   input_file.close();
      
   //Perform the double probe curve fit using non-linear least squares
   std::cout << "Performing curve fit..." << std::endl;
   if(IVFit2NLLS(Ii, Vi, Max, Tol, FitParams)){
      
      std::cout << "Curve fit successful!" << std::endl;
      
   }else{
    
      std::cout << "Curve fit failed" << std::endl;
      
   }
   
   //Print the fitted parameters
   std::cout.precision(3);
   std::cout << "Final fit parameters: " << std::endl;
   std::cout << " Ion saturation current [A]  : " << FitParams.Isat << std::endl;
   std::cout << " Electron temperature   [eV] : " << FitParams.Te << std::endl;
   
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
      double col1 = 0.0,
             col2 = 0.0;
     
      for(std::vector<double>::iterator itr = Vi.begin(); itr != Vi.end();
                                                                   ++itr){
            
            col1 = *itr;
            col2 = Iv(col1, FitParams.Isat, FitParams.Te);
            output_file << col1 << " ";
            output_file << col2 << std::endl;
         
      }
     
   }else{
     
      std::cerr << "Error opening file:" << output_filename_s.c_str();
      std::cerr << std::endl;
      return (-1);
      
   }//Done writing output file
   
   output_file.close();
      
 
std::cout << "-- END DoubleProbeAnalysis --" << std::endl;
return(0);

}

/************************************************************************/
void print_usage(){
   
   std::cout << "Usage:" << std::endl;
   std::cout << "bin/DoubleProveAnalysis -f <filename>" << std::endl;
   std::cout << "bin/DoubleProbeAnalysis -f ExampleData/ExampleData.dat";
   std::cout << std::endl;
   
}
