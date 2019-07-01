README for IFOS3D
=================

IFOS3D (Inversion of Full observed seismograms) 
is a 3D elastic full waveform inversion code. The inversion problem is solved 
by a conjugate gradient method. The gradients are computed in the frequency domain 
for discrete frequencies with the adjoint method. The forward modeling is done 
in the time domain with a finite difference code on a standard staggered grid.


Installation instructions can be found in the file INSTALL.


A detailed documentation is provided in manual_IFOS3D.pdf.


To get started with this code try to run the toy example that is described 
in the documentation (manual_IFOS3D.pdf).


If you use this code for your own research please cite at least one article
written by the developers of the package, e.g.
Butzer, S., Kurzmann, A., & Bohlen, T., 2013. 3D elastic full-waveform inversion of small-
scale heterogeneities in transmission geometry, Geophysical Prospecting, 61(6), 1238-1251.


Contact:
see file AUTHORS
