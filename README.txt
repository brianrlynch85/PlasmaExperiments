# -------------------------------------------------------------------------------------
#
#                            README for DoubleProbeAnalysis
#                                        V 0.01
#
#                            (c) Brian Lynch February, 2015
#
# -------------------------------------------------------------------------------------

This software was written for educational purposes. It is intended to be a walkthrough
of performing a simple nonlinear least squares curve fit on double Langmuir probe data
taken from a plasma physics experiment. For more information about nonlinear least
squares and the Langmuir probe diagnostic, see the following links:

http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html
http://en.wikipedia.org/wiki/Langmuir_probe

In this code, LAPACK libraries are used with the -llapack link in the Makefile.
For more information about the LAPACK linear algebra library, see the following link:

   http://www.netlib.org/lapack/
   
A copy of the LAPACK modified BSD license is contained in the subfolder
titled 'matrix_utils'. You will likely need to execute the command
"sudo apt-get install liblapack-dev" to install the liblapack libraries.

This software works on Ubuntu 12.04 & 14.04 and has not been verified on other
operating systems.
   possible make options are:
      "make clean"            to clean up all .o and temp files
      "make mrclean"          to clean up all .o, temp files, and executables
      "make all"              to compile everything
      
Example calling commands (using example data provided in the ExampleData folder):
   bin/DoubleProbeAnalysis -f <inputfilename>
   bin/DoubleProbeAnalysis -f ExampleData/ExampleData.dat
   
When using the example data, you should get the following terminal output:

   > bin/DoubleProbeAnalysis -f ExampleData/ExampleData.dat
   > -- BEGIN DoubleProbeAnalysis --
   > Input Filename: ExampleData/ExampleData.dat
   > Initial fit parameters: 
   >  Ion saturation current [A]  : 3.3e-06
   >  Electron temperature   [eV] : 3
   > Reading IV data...
   > Performing curve fit...
   >  R^2         : 5.81e-09
   >  # iterations: 8
   > Curve fit successful!
   > Final fit parameters: 
   >  Ion saturation current [A]  : 8.54e-06
   >  Electron temperature   [eV] : 21.1
   > Writing fit data to file: ExampleData/ExampleData_fit.dat
   > -- END DoubleProbeAnalysis --

Inside the folder ExampleData, I included a png file titled "ExampleDataPlot.png"
as an example of the output produced by the code. It should be noted that the fit
is not perfect. In particular, the electron temperature is a bit high (in reality
it is more like 4eV). The descrepancy between the data and curve fit is due to the
idealization made when deriving the I(V) fit. Better, slightly more complicated,
models include effects geometric and plasma sheath effects.