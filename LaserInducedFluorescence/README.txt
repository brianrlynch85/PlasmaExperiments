README for LaserInducedFluorescence
===================================
(c) Brian Lynch March, 2015
---------------------------

This software was written for educational purposes. It is intended to be a
walkthrough of performing a simple nonlinear least squares curve fit on
Laser Induced Flurescence data taken from a plasma physics experiment.
For more information about nonlinear least squares and the Laser Induced
Fluorescence diagnostic, see the following links:

http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html
http://en.wikipedia.org/wiki/Laser-induced_fluorescence

In this code, LAPACK libraries are used with the -llapack link in the
Makefile. For more information about the LAPACK linear algebra library, see
the following link:

   http://www.netlib.org/lapack/
   
A copy of the LAPACK modified BSD license is contained in the subfolder
titled 'matrix_utils'. You will likely need to execute the command
"sudo apt-get install liblapack-dev" to install the liblapack libraries.

You may also need to install cmake. To install cmake, execute the command
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
         build/bin/LIFAnalysis -f <inputfilename>
         build/bin/LIFAnalysis -f ExampleData/ExampleData.dat
   
When using the example data, you should get the following terminal output:

> bin/LIFAnalysis -f ../ExampleData/ExampleData.dat
> -- BEGIN lif_analysis --
> Input Filename: ../ExampleData/ExampleData.dat
> Initial fit parameters: 
>  Rest Wavelength        [nm]  : 668.6138
>  Sigma^2               [nm^2] : 6e-07
>  Amplitude               []   : 4
>  Background              []   : 0.5
> Reading data...
> Performing curve fit...
>  R^2         : 1.409663e-09
>  # iterations: 6
> Curve fit successful!
> Final fit parameters: 
>  Rest Wavelength        [nm]  : 668.6137
>  Sigma^2               [nm^2] : 5.207346e-07
>  Amplitude               []   : 3.740877
>  Background              []   : 0.4728978
> Writing fit data to file: ../ExampleData/ExampleData_fit.dat
> -- END lif_analysis --

Inside the folder ExampleData, you will find a .png file titled
"ExampleDataPlot.png". It is an example of the output produced by the
code. Note that the fit is not perfect. The known transition is at
a wavelength of 668.6138 nm. From the width of the distribution, we
could obtain the ion temperature. From the shift of the mean of the
distribution, we could obtain the ion drift velocity.