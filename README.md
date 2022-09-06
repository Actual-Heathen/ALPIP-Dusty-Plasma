# ALPIP-Dusty-Plasma

Dust particles have a ubiquitous presence in space, being found in interplanetary and interstellar space. Dust is found in nebulae which are giant clouds of dust and plasma. Dust particles in the plasma are able to couple with the plasma. This coupling of dust modifies the linear dispersion relations of waves within the plasma. Furthermore, charged dust particles are subjected to the Lorentz force from the ambient magnetic field. In addition, they experience self-gravitation due to the mass of the particles within them. This gravitation leads to instabilities such as Jeans instability. Jeans instability is caused when the dust pressure is unable to counteract the gravitational forces and the particles start to self-contract. In the code, we show the linear to non-linear evolution of the dust particles over time with a background magnetic field using the Particle-In-Cell (PIC) simulation method. The gravitational potential is solved by the Poisson equation. The interaction with the background magnetic field is taken into account in the calculation of the Lorentz force with a python frontend for visualization.

# Dependencies

To compile and run the code the following should be installed:
The C++ code requires the [FFTW library](https://www.fftw.org/) to perform the necessary fast fourier transforms.
The python script requires matplotlib and pandas to be installed.

# Acknowledgments

This work is supported by the NSF EPSCoR RII-Track-1 Cooperative Agreement OIA-1655280

