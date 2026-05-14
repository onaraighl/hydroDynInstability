# TurbulenceModel

This directory contains a simple code to compute the mean velocity profile in a channel in case of a turbulent channel flow.  The code uses the turbulence model for channel flows describedin <b>Chapter 11</b> of the reference textbook.

The main code is ``get_U_singlephase.m``  The single input is the Renolds number $ReG$.  The outputs are:

* zG - z-coordinate
* uG - average $u$-velcocity, as a function of $z$.

These can be visualized, e.g.

``plot(zG,uG)``

Further outputs:

* Re_star - Reynolds number based on friction velocity
* ReG - redundant, as this is an input
* tauG - Reynolds stress, again $z$-dependent
* kG - turbulent kinetic energy, again $z$-depenent

Sample results using this code are shown below.

![Cartoon](cartoon1.png1)
