# ------------------------------------------------------------------------
#
#                            README for DoubleProbeAnalysis
#                                        V 0.01
#
#                            (c) Brian Lynch February, 2015
#
# ------------------------------------------------------------------------

This software was written for educational purposes. It is intended to be a
walkthrough of performing a simple nonlinear least squares curve fit on
double Langmuir probe data taken from a plasma physics experiment. For
more information about nonlinear least squares and the Langmuir probe
diagnostic, see the following links:

http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html
http://en.wikipedia.org/wiki/Langmuir_probe

In this code, LAPACK libraries are used with the -llapack link in the
Makefile. For more information about the LAPACK linear algebra library, see
the following link:

   http://www.netlib.org/lapack/
   
A copy of the LAPACK modified BSD license is contained in the subfolder
titled 'matrix_utils'. You will likely need to execute the command
"sudo apt-get install liblapack-dev" to install the liblapack libraries.

You may also need to install cmake depending whether or not you choose
to use the older file Makefile.Old. To install cmake, execute the command
"sudo apt-get install cmake".

This software works on Ubuntu 12.04 & 14.04 and has not been verified on other
operating systems.

   possible cmake options are (will put the executables in build/bin):
      "mkdir build"
      "cd build"
      "cmake ../"
      "make"
      "cd ../"
      Now you have done and out of source build, which leaves the original
      source directories clean.
      
      Example calling commands (using example data in the ExampleData folder):
         build/bin/DoubleProbeAnalysis -f <inputfilename>
         build/bin/DoubleProbeAnalysis -f ExampleData/ExampleData.dat
   

   possible make options are (will put the executables in bin):
      "make -f Makefile.Old clean"            to clean up all .o and temp files
      "make -f Makefile.Old mrclean"          to clean up all .o, temp files, and executables
      "make -f Makefile.Old all"              to compile everything
      
      Example calling commands (using example data in the ExampleData folder):
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

Inside the folder ExampleData, you will find a .png file titled
"ExampleDataPlot.png". It is an example of the output produced by the
code. Note that the fit is not perfect. In particular, the electron
temperature is a bit high (in reality it is more like 8eV).
The descrepancy between the data and curve fit is due to the idealizations
made when deriving the simple 2 parameter I(V) fit. There exists other,
slightly more complicated, models which include geometric and plasma
sheath effects.