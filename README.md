# What is IFOS3D?

**IFOS3D** (**I**nversion of **F**ull **O**bserved **S**eismograms) is a 3D elastic full waveform inversion code.  
The inversion problem is solved by a conjugate gradient method and the gradients are computed in the frequency domain with the adjoint method.  
The forward modeling is done by a time domain finite difference scheme. 

This branch is aimed to accelerate the inversion process by first fitting the low frequency samples on coarse grids, and gradually refine the grid towards higher frequencies.  

IFOS3D is the reverse (inverse) of our 3D finite difference forward solver [**SOFI3D**](https://git.scc.kit.edu/GPIAG-Software/SOFI3D).
