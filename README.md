Plasma Experiments
==================
(c) Brian Lynch June, 2016

This repo contains a collection of plasma physics data and its corresponding
analysis code. Each project is self contained with its own README.md file.

To-do:
######

* I should probably combine the project more efficiently so that they use
the same matrix_utils instead of them each having their own version of the
exact same matrix utilities.

* The nonlinear least squares methods are essentially the same (just with
different Jacobians). I should probably abstract these routines a bit and
use a more object oriented/C++ approach to avoid redundancy.

* Latex/doxygen documentation of the code and tutorials on the Physics
contained in the data.
